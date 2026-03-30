function p = getSimParams()
% =========================================================
% Function: getSimParams.m
% Purpose : Single source of truth for Table I simulation
%           parameters. All four entry-point scripts call
%           this function and override only what they change.
% Paper   : Sensors-100615-2025
%           "Covariance-Guided DFT Beam Selection for
%            Beamspace ESPRIT in Hybrid mmWave Sensor Arrays"
% Author  : Rifat Volkan Senyuva
% =========================================================
%
% USAGE
%   p = getSimParams();
%   % Use p.M, p.N, p.NRF, p.mu_true, etc.
%   % Override script-specific fields after the call, e.g.:
%   %   p.ASNR = -4:1:6;    % A1 uses a narrower ASNR range
%
% OUTPUT FIELDS
%   Array / hardware
%     p.M          : number of array elements (32)
%     p.N          : number of snapshots per trial (100)
%     p.NRF        : number of RF chains / virtual subarray size (12)
%     p.d          : number of sources (3)
%
%   Source model (Table I)
%     p.pows       : source powers [0.95 0.50 0.10]
%     p.mu_true    : true spatial frequencies [rad] (3x1)
%     p.R_source   : diagonal source covariance (dxd)
%     p.Psqrt      : sqrt(pows/2) diagonal — for complex Gaussian draws
%
%   Beam selection
%     p.W          : DFT sector window width in beams (4)
%     p.zeta       : sectorize_dft_beams offset parameter (1)
%     p.Kg         : beam budget per sector for fine stage (2)
%     p.alpha      : conditioning penalty weight in score (0.5)
%     p.eta        : steering tolerance in score (0, inactive for unitary DFT)
%     p.k_thr      : failure-rate threshold multiplier (3)
%
%   Simulation control
%     p.ASNR       : default ASNR sweep [-10:1:20] dB  (Main, S1 override)
%     p.totalIter  : Monte Carlo trials per condition (1e4)
%     p.rng_seed   : Threefry seed for RandStream (5489)
%
%   Selector options struct (matches select_*_from_sectorized_cov)
%     p.opts_sel   : struct with eig_tol_scale, tikhonov_gamma,
%                    alpha, normalize_score, prune_topq
%
%   Derived quantities (precomputed here for convenience)
%     p.Atrue      : M x d true steering matrix
%     p.R_signal   : M x M true signal covariance
%     p.phase_k    : M x 1 FFT phase correction vector
%     p.maskIdx    : 1 x NRF centered-contiguous subarray indices
%     p.mIdxFlip   : NRF:-1:1  (FBA reflection index)
%
% NOTES
%   * paperCodeA2 (sector-edge test) uses d=2, pows=[0.95 0.50],
%     mu_base=[-2.1;0.5], W=6, ASNR_set=[3,6]. These deviations
%     are set locally in paperCodeA2.m after calling getSimParams.
%   * Timing benchmarks in data/A1_timing_table.csv were generated
%     on an older desktop. Timings on the current machine will be
%     faster; compare RMSE_rad only, not wall-clock columns.
%   * The filename estimate_powers_course.m contains a typo
%     ('course' for 'coarse'). It is preserved as-is for compatibility
%     with existing results; rename planned for Phase 2 refactoring.
% =========================================================

%% ---- Array / hardware --------------------------------------------------
p.M   = 32;
p.N   = 100;
p.NRF = 12;

%% ---- Source model (Table I) -------------------------------------------
p.pows    = [0.95, 0.50, 0.10];
p.mu_true = [-2.1; 0.5; 2.5];          % spatial frequencies [rad]
p.d       = numel(p.mu_true);

p.R_source = diag(p.pows);
p.Psqrt    = diag(sqrt(p.pows / 2));   % complex Gaussian amplitude scaling

%% ---- Beam selection parameters ----------------------------------------
p.W     = 4;          % sector window width [beams]  — A2 overrides to 6
p.zeta  = 1;          % sectorize_dft_beams offset parameter
p.Kg    = 2;          % beam budget per sector (K_FINE in A1, kf_list in S1)
p.alpha = 0.5;        % conditioning penalty weight
p.eta   = 0;          % steering tolerance (0 => inactive for unitary DFT)
p.k_thr = 3;          % failure-rate threshold multiplier (kThr in paper)

%% ---- Simulation control -----------------------------------------------
p.ASNR      = -10:1:20;   % default sweep (A1 narrows to -4:6; S1 uses [0 5 10 15])
p.totalIter = 1e4;
p.rng_seed  = 5489;        % Threefry seed — do NOT change; fixes Monte Carlo stream

%% ---- Default selector options struct ----------------------------------
% Used by select_adjacent_pairs_from_sectorized_cov and
% select_beamsets_from_sectorized_cov (paperCodeMain uses defaults,
% A1/S1 override prune_topq).
p.opts_sel = struct( ...
    'eig_tol_scale',    1e-6, ...
    'tikhonov_gamma',   0,    ...
    'alpha',            0.5,  ...
    'normalize_score',  true, ...
    'prune_topq',       0     );   % 0 = off; A1 overrides to 4

%% ---- Derived quantities (computed once here) --------------------------
p.Atrue   = exp(1i * (0:p.M-1).' * p.mu_true.');   % M x d steering matrix
p.R_signal = p.Atrue * p.R_source * p.Atrue';      % M x M signal covariance

% FFT phase correction (applied row-wise after fft(...,M,1))
p.phase_k  = exp(1i * (p.M - 1) * pi * (0:p.M-1).' / p.M);  % M x 1

% Virtual subarray indices (requires centered_contiguous_mask on path)
p.maskIdx  = centered_contiguous_mask(p.M, p.NRF);   % 1 x NRF
p.mIdxFlip = p.NRF:-1:1;                             % FBA reflection index

end
