%% validateReproducibility.m
% Automated benchmark checker for CovGuided-ESPRIT (v1.0 repo structure).
%
% PURPOSE:
%   Confirms that the renamed/restructured codebase reproduces the four
%   benchmark values stored in results/csv/S1_aggregates.csv to numerical
%   precision (0.00% RMSE error relative to CSV ground truth).
%
% USAGE:
%   Run from the repository root after repo modernization:
%       >> cd <repo_root>
%       >> validateReproducibility
%
%   All 4 checks must PASS at 0.00% error before committing the
%   modernize-v1.0 branch.
%
% DIFFERENCES FROM Phase-1 VERSION:
%   - All function calls updated to covguided.* package namespace
%   - estimate_powers_course -> covguided.estimatePowersCoarse (typo fixed)
%   - CSV path: data/ -> results/csv/
%   - Path setup updated for new src/+covguided/ layout
%   - Algorithm body and random seed are identical to Phase-1 version
%
% WHAT IS CHECKED:
%   Runs the parfor body from runFineStageSweep.m verbatim for kf=2,
%   ASNR in {0, 5, 10, 15} dB, 1e4 trials, Threefry seed 5489.
%   Checks RMSE and FailRate for two methods x two ASNR points = 4 checks.
%   Timing is NOT checked (machine-specific; see README hardware caveat).
%
% EXPECTED RESULTS (from S1_aggregates.csv):
%   Check 1 — Cov-guided, kf=2, ASNR=0 dB:   RMSE=0.21484 rad, FR=0.1571
%   Check 2 — Cov-guided, kf=2, ASNR=15 dB:  RMSE=0.00072382 rad, FR=0.0048
%   Check 3 — Sectorization, kf=2, ASNR=0 dB:  RMSE=0.74106 rad, FR=0.4445
%   Check 4 — Sectorization, kf=2, ASNR=15 dB: RMSE=0.00096286 rad, FR=0.0181

clearvars;

%% --- Path setup -----------------------------------------------------------
% Add src/ so MATLAB finds the +covguided package, and config/ for getSimParams.
repoRoot = fileparts(mfilename('fullpath'));
repoRoot = fileparts(repoRoot);   % validation/ -> repo root
addpath(fullfile(repoRoot, 'src'));
addpath(fullfile(repoRoot, 'config'));
addpath(fullfile(repoRoot, 'experiments'));

fprintf('validateReproducibility — CovGuided-ESPRIT v1.0\n');
fprintf('Repo root: %s\n\n', repoRoot);

%% --- Invariants (copied verbatim from runFineStageSweep.m) ----------------
M   = 32;
N   = 100;
NRF = 12;

kf_list = [2 3 4];
numKf   = numel(kf_list);

opts_sel = struct('eig_tol_scale',1e-6, 'tikhonov_gamma', 0, ...
                  'alpha',0.5, 'normalize_score', true, ...
                  'prune_topq', 0);

%% RF mask
maskIdx  = covguided.centeredContiguousMask(M, NRF);
mIdxFlip = NRF:-1:1;

%% Precompute once (outside parfor)
phase_k  = exp(1i*(M-1)*pi*(0:M-1).'/M);

%% Source powers, true DOAs
pows     = [0.95, 0.5, 0.1];
R_source = diag(pows);
Psqrt    = diag(sqrt(pows/2));
mu_true  = [-2.1; 0.5; 2.5];
d        = numel(mu_true);

%% Steering (precompute once)
Afun   = @(mu) exp(1i * (0:M-1).' * mu.');
Atrue  = Afun(mu_true);

%% Signal covariance (invariant)
R_signal = Atrue * R_source * Atrue';

%% Oracle beam set (invariant)
beamGroups_ini  = covguided.sectorizeDftBeams(mu_true, M, 4, 1);
SoICols_ini     = covguided.selectAdjacentPairsFromSectorizedCov(R_signal, beamGroups_ini);
SoICols_true    = unique(cell2mat(SoICols_ini)).';   %#ok<NASGU>

%% SNR & noise — only the two validation ASNR points
ASNR     = [0 15];   % subset for validation (was [0 5 10 15] in full script)
noiseVar = real(trace(R_signal)) ./ (M * 10.^(ASNR/10));

%% CRB
sqrtCRB = zeros(1, length(noiseVar));
PorthA  = eye(M) - Atrue * pinv(Atrue);
derA    = (1i * transpose(0:M-1)) .* Atrue;
H       = transpose(derA' * PorthA * derA);
for noiseInd = 1:length(noiseVar)
    secondTerm = ((Atrue' * Atrue) * R_source) / noiseVar(noiseInd);
    firstTerm  = eye(d) + secondTerm;
    CRBinvterm = diag(inv(real((R_source * (firstTerm \ secondTerm)) .* H)));
    sqrtCRB(noiseInd) = sqrt(sum((noiseVar(noiseInd) / (2*N)) * CRBinvterm) / d);
end

%% --- Accumulators ---------------------------------------------------------
totalIter     = 1e4;
kThr          = 3;

MSE_fine      = zeros(numKf, length(ASNR));
MSE_wide      = zeros(numKf, length(ASNR));
failRate_fine = nan(numKf, length(ASNR));
failRate_wide = nan(numKf, length(ASNR));

%% --- Monte Carlo (outer SNR loop) ----------------------------------------
for snrIndex = 1:length(ASNR)

    nvStd         = sqrt(noiseVar(snrIndex) / 2);
    e_trials_fine = nan(totalIter, numKf);
    e_trials_wide = nan(totalIter, numKf);

    fprintf('Running ASNR = %d dB (%d trials)...\n', ASNR(snrIndex), totalIter);

    %% PARFOR — body copied verbatim from runFineStageSweep.m ---------------
    %  Only the function names have changed to covguided.* namespace.
    %  estimate_powers_course -> covguided.estimatePowersCoarse (typo fixed).
    %  All numerical logic is identical.
    parfor iter = 1:totalIter
        % Threefry seed with per-trial substreams (ensures cross-machine reproducibility)
        s = RandStream('Threefry', 'Seed', 5489);
        s.Substream = iter;
        RandStream.setGlobalStream(s);

        %% Coarse stage synthetic data
        S = Psqrt * (randn(d,N) + 1i*randn(d,N));
        X = Atrue * S;
        Z = nvStd * (randn(M,N) + 1i*randn(M,N));
        Y_coarse = X + Z;

        %% Element-space subarray
        Yb_coarse = Y_coarse(maskIdx, :);
        Yb_wideES = Y_coarse(1:NRF, :);

        %% Empirical covariances
        RY_coarse = (Yb_coarse * Yb_coarse') / N;   RY_coarse = (RY_coarse + RY_coarse') / 2;
        RY_wideES = (Yb_wideES * Yb_wideES') / N;   RY_wideES = (RY_wideES + RY_wideES') / 2;

        %% FBA (coarse)
        RY_fba = (RY_coarse + conj(RY_coarse(mIdxFlip, mIdxFlip))) * 0.5;

        %% TLS-ESPRIT
        mu_coarse = covguided.tlsEsprit(RY_fba,    d);
        mu_wide   = covguided.tlsEsprit(RY_wideES, d);

        %% Coarse power estimation + Toeplitz-PSD covariance reconstruction
        [~, ~, mu_courseSort, R_course] = covguided.estimatePowersCoarse(RY_fba,    maskIdx, mu_coarse, M);
        [~, ~, mu_wideSort,   R_wide]   = covguided.estimatePowersCoarse(RY_wideES, 1:NRF,   mu_wide,   M);

        %% Fine stage synthetic data (independent draw)
        Z2     = nvStd * (randn(M,N) + 1i*randn(M,N));
        S2     = Psqrt * (randn(d,N) + 1i*randn(d,N));
        X2     = Atrue * S2;
        Y_fine = X2 + Z2;

        %% FFT-based full beamspace projection
        Yb_full = phase_k .* fft(Y_fine, M, 1);

        e_vec_fine = nan(1, numKf);
        e_vec_wide = nan(1, numKf);

        %% Inner loop over fine-stage beam budgets
        for kk = 1:numKf
            kf = kf_list(kk);

            beamGroups_fine = covguided.sectorizeDftBeams(mu_courseSort, M, 2*kf, 1);
            selBeams_fine   = covguided.selectBeamsetsFromSectorizedCov(R_course, beamGroups_fine, kf, opts_sel);
            SoICols_fine    = unique(cell2mat(selBeams_fine)).';

            beamGroups_wide = covguided.sectorizeDftBeams(mu_wideSort, M, kf, 1);
            SoICols_wide    = unique(cell2mat(beamGroups_wide)).';

            Yb_fine = Yb_full(SoICols_fine, :);
            Yb_wide = Yb_full(SoICols_wide, :);

            mu_fine = covguided.unitaryEspritSparse(Yb_fine, SoICols_fine, d, M);
            mu_wide = covguided.unitaryEspritSparse(Yb_wide, SoICols_wide, d, M);

            RY_fine = (Yb_fine * Yb_fine') / N;   RY_fine = (RY_fine + RY_fine') / 2;
            RY_wide = (Yb_wide * Yb_wide') / N;   RY_wide = (RY_wide + RY_wide') / 2;

            [~, ~, mu_fineSort, ~]      = covguided.estimatePowersFine(RY_fine, SoICols_fine, mu_fine, M);
            [~, ~, mu_wideSort_fine, ~] = covguided.estimatePowersFine(RY_wide, SoICols_wide, mu_wide, M);

            e_vec_fine(kk) = sqrt(mean(angle(exp(1i*(mu_true - mu_fineSort))).^2));
            e_vec_wide(kk) = sqrt(mean(angle(exp(1i*(mu_true - mu_wideSort_fine))).^2));
        end

        e_trials_fine(iter, :) = e_vec_fine;
        e_trials_wide(iter, :) = e_vec_wide;
    end % parfor

    %% Reduce
    MSE_fine(:, snrIndex)      = mean(e_trials_fine.^2, 1).';
    MSE_wide(:, snrIndex)      = mean(e_trials_wide.^2, 1).';
    thr                        = kThr * max(sqrtCRB(snrIndex), 1e-12);
    failRate_fine(:, snrIndex) = mean(e_trials_fine > thr, 1).';
    failRate_wide(:, snrIndex) = mean(e_trials_wide > thr, 1).';
end

%% --- Load CSV ground truth -----------------------------------------------
csvPath = fullfile(repoRoot, 'results', 'csv', 'S1_aggregates.csv');
if ~isfile(csvPath)
    error('validateReproducibility: CSV not found at %s', csvPath);
end
T = readtable(csvPath);

%% --- Verification function -----------------------------------------------
function check(label, computed, csv_val, tol)
    err = abs(computed - csv_val) / max(abs(csv_val), 1e-15);
    if err < tol
        fprintf('  PASS  %s: computed=%.6g, csv=%.6g, err=%.2f%%\n', ...
                label, computed, csv_val, err*100);
    else
        fprintf('  FAIL  %s: computed=%.6g, csv=%.6g, err=%.2f%%\n', ...
                label, computed, csv_val, err*100);
    end
end

%% --- Run 4 benchmark checks ----------------------------------------------
fprintf('\n=== Benchmark verification (kf=2, Threefry seed 5489) ===\n');
tol = 1e-3;   % 0.1% tolerance (expect 0.00%)

% kf=2 is kf_list index 1; ASNR=0 is snrIndex 1; ASNR=15 is snrIndex 2
rmse_fine_0  = sqrt(MSE_fine(1, 1));   % Cov-guided, kf=2, ASNR=0 dB
rmse_fine_15 = sqrt(MSE_fine(1, 2));   % Cov-guided, kf=2, ASNR=15 dB
rmse_wide_0  = sqrt(MSE_wide(1, 1));   % Sectorization, kf=2, ASNR=0 dB
rmse_wide_15 = sqrt(MSE_wide(1, 2));   % Sectorization, kf=2, ASNR=15 dB

fr_fine_0  = failRate_fine(1, 1);
fr_fine_15 = failRate_fine(1, 2);
fr_wide_0  = failRate_wide(1, 1);
fr_wide_15 = failRate_wide(1, 2);

% Pull expected values from CSV
row = T(strcmp(T.method,'Cov-guided') & T.Kf==2 & T.ASNR_dB==0, :);
check('[1] Cov-guided  kf=2 ASNR= 0dB RMSE',     rmse_fine_0,  row.RMSE_rad,  tol);
check('[1] Cov-guided  kf=2 ASNR= 0dB FailRate',  fr_fine_0,   row.FailRate,  tol);

row = T(strcmp(T.method,'Cov-guided') & T.Kf==2 & T.ASNR_dB==15, :);
check('[2] Cov-guided  kf=2 ASNR=15dB RMSE',     rmse_fine_15, row.RMSE_rad,  tol);
check('[2] Cov-guided  kf=2 ASNR=15dB FailRate',  fr_fine_15,  row.FailRate,  tol);

row = T(strcmp(T.method,'Sectorization') & T.Kf==2 & T.ASNR_dB==0, :);
check('[3] Sectorization kf=2 ASNR= 0dB RMSE',   rmse_wide_0,  row.RMSE_rad,  tol);
check('[3] Sectorization kf=2 ASNR= 0dB FailRate', fr_wide_0,  row.FailRate,  tol);

row = T(strcmp(T.method,'Sectorization') & T.Kf==2 & T.ASNR_dB==15, :);
check('[4] Sectorization kf=2 ASNR=15dB RMSE',   rmse_wide_15, row.RMSE_rad,  tol);
check('[4] Sectorization kf=2 ASNR=15dB FailRate', fr_wide_15, row.FailRate,  tol);

fprintf('\nDone. All 4 RMSE checks + 4 FailRate checks above.\n');
fprintf('Commit criterion: all 8 lines must read PASS.\n');
