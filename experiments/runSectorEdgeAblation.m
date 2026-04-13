clearvars;
addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'config'));

%% --- Output directory setup (repo-relative) -------------------
repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
figDir   = fullfile(repoRoot, 'results', 'figures');
csvDir   = fullfile(repoRoot, 'results', 'csv');
if ~exist(figDir,'dir'), mkdir(figDir); end
if ~exist(csvDir,'dir'), mkdir(csvDir); end
%% ---------------------------------------------------------------

% ==== A2: Sector-edge stress test =========================================
M          = 32;
N          = 100;
NRF_coarse = 12;                   % coarse=12 (ES mask)
NRF_fine   = 6;                    % fine=6 beams

W          = 6;                    % sector window width in beams
K          = 17;                   % 15–21 grid points across the sector width

maskIdx    = covguided.centeredContiguousMask(M, NRF_coarse);       % reuse your mask
mIdxFlip   = NRF_coarse:-1:1;                               % for FBA reflect (coarse size)

%% Steering (precompute once)
Afun    = @(mu) exp(1i * (0:M-1).' * mu.'); 
phase_k = exp(1i * (M-1) * pi * (0:M-1).'/M );      % size Mx1, for row-wise scaling

%% Choose a boundary between beams k and k+1 to stress
%k_bdy   = 8;                                 % representative interior boundary
%gamma   = @(b) (2*pi/M)*b - pi;              % beam center ([-pi,pi))
%mu_edge = 0.5*(gamma(k_bdy) + gamma(k_bdy+1));
%w       = W*(2*pi/M);
%delta   = linspace(-w/2, +w/2, K);           % sweep across the sector

% Beam and sector geometry
beam_spacing = 2 * pi / M;                                   % between adjacent DFT beams in mu
w            = W * beam_spacing;                          % sector width in mu, still defined

% Choose a boundary between beams k_bdy and k_bdy+1 to stress (near broadside)
gamma       = @(b) (2*pi/M)*b - pi;                      % beam center ([-pi,pi))
k_bdy       = M/2 - 1;                                   % e.g., 15 for M=32
mu_edge     = 0.5 * (gamma(k_bdy) + gamma(k_bdy+1));

% Sweep only within ±alpha beam spacings around the edge
alpha        = 1.5;                                      % in beam spacings; try 1.0–1.5
delta        = linspace(-alpha * beam_spacing, +alpha * beam_spacing, K);


%% Freeze the other sources as in your defaults; move exactly one
%pows     = [0.95, 0.5, 0.1];
pows     = [0.95, 0.5];                        % 2 paths now
R_source = diag(pows);
sumPows  = sum(pows);
d        = numel(pows);
Psqrt    = diag(sqrt(pows/2));
%mu_base  = [-2.1; 0.5; 2.5];                 % your default
mu_base  = [-2.1; 0.5];                       % 2 x 1
idx_move = 2;                                 % move this DOA to the boundary

ASNR_set = [3, 6];                            % dB
noiseVar = sumPows ./ 10.^(ASNR_set/10);

%% number of methods
totalIter  = 1e4;
numMethods = 2;                               % 1=cov-guided, 2=sectorization

%% Accumulators: failure rate (%) vs delta, per SNR and method
Fail_vs_delta = zeros(numel(ASNR_set), numMethods, K);
MSE           = zeros(numel(ASNR_set), numMethods, K);

% Failure threshold factor
kThr          = 3;
thrVal        = zeros(numel(ASNR_set), K);
sqrt_CRB       = zeros(numel(ASNR_set), K);



%% ========================= MAIN LOOPS =========================
fprintf('runSectorEdgeAblation: %d SNR points x %d delta points = %d ticks total\n', ...
        numel(ASNR_set), K, numel(ASNR_set)*K);
totalTicks = numel(ASNR_set) * K;
ticksDone  = 0;
tStart     = tic;

for sIdx = 1:numel(ASNR_set)

    nv    = noiseVar(sIdx);
    nvStd = sqrt(nv/2);

    % --- run the Monte Carlo engine per delta point -------------
    for t = 1:K

        % Per-δ local accumulators (parfor-friendly)
        e_trials          = zeros(totalIter, numMethods);
        fail_flags_fine   = false(totalIter, 1);
        fail_flags_sector = false(totalIter, 1);

        %% True geometry at this δ
        mu_true           = mu_base;
        mu_true(idx_move) = wrapToPi(mu_edge + delta(t));                     % place at boundary offset
        Atrue             = Afun(mu_true);

        %% --- Per-source CRB in mu (spatial frequency), uncorrelated Gaussian signals ---
        PorthA            = eye(M) - Atrue*pinv(Atrue);                       % noise subspace projector
        D                 = (1i*(0:M-1).').*Atrue;                            % Mxd, dA/dmu
        G                 = real(D' * PorthA * D);                            % dxd, geometry term (Hermitian)
        Wsnr              = R_source / nv;                                    % diag of per-source SNRs (linear)
        Fmu               = 2 * N * (G .* Wsnr);                              % Hadamard weighting for uncorrelated sources
        CRBmu             = inv(Fmu);                                         % dxd
        sqrt_CRB(sIdx,t)  = sqrt(mean(diag(CRBmu)));                          % sqrt(mean variance across sources)
        %
        sigma_mu_src      = sqrt(max(diag(CRBmu), 0));                        % 1xd, nonnegative
        thr_src           = transpose(kThr * max(sigma_mu_src, 1e-12));   
        thrVal(sIdx, t)   = max(thr_src);
        

        %% PARFOR over Monte-Carlo iterations
        parfor iter = 1:totalIter
            % Stream per iteration (deterministic & independent)
            s = RandStream('Threefry', 'Seed', 5489);
            s.Substream = iter;
            RandStream.setGlobalStream(s);

            % --- Coarse stage synthetic data (element space) ---
            S        = Psqrt * (randn(d,N) + 1i*randn(d,N));
            X        = Atrue * S;

            % noise and received
            Z        = nvStd * (randn(M,N) + 1i*randn(M,N));
            Y_coarse = X + Z;

            % element-space subarray "beamforming"
            Yb_coarse = Y_coarse(maskIdx, :);          % NRF_coarse x N (contiguous mask)

            % empirical covariances (Hermitian symmetrization)
            RY_coarse = (Yb_coarse * Yb_coarse') / N;   RY_coarse = (RY_coarse + RY_coarse') / 2;

            % FBA (coarse)
            RY_fba    = 0.5 * (RY_coarse + conj(RY_coarse(mIdxFlip, mIdxFlip)));

            % TLS-ESPRIT on coarse/wideES covariances
            mu_coarse = covguided.tlsEsprit(RY_fba,     d);

            % Coarse power estimates (unchanged)
            [~, ~, mu_courseSort, R_course] = covguided.estimatePowersCoarse(RY_fba,     maskIdx,      mu_coarse, M);
            
            % Single sectorization call with W:
            beamGroups    = covguided.sectorizeDftBeams(mu_courseSort, M, W, 1);

            % Cov-guided: pick adjacent pair(s) from the *same* coarse covariance
            selPairs      = covguided.selectAdjacentPairsFromSectorizedCov(R_course, beamGroups);
            SoICols_fine  = pad_truncate_to_nrf(unique(cell2mat(selPairs)).', NRF_fine, M);
            
            % Sectorization baseline: fixed window from the same beamGroups
            SoICols_wide  = pad_truncate_to_nrf(unique(cell2mat(beamGroups)).', NRF_fine, M);

            % --- Fine stage synthetic data (element space) ---
            Z2     = nvStd * (randn(M,N) + 1i*randn(M,N));
            S2     = Psqrt * (randn(d,N) + 1i*randn(d,N));
            X2     = Atrue * S2;
            Y_fine = X2 + Z2;

            % --- FFT path to DFT beamspace (big win) ---
            Yb_full = phase_k .* fft(Y_fine, M, 1);

            % Build Yb_* by row-indexing with SoI
            Yb_fine = Yb_full(SoICols_fine, :);
            Yb_wide = Yb_full(SoICols_wide, :);

            % Run the same fine-stage ESPRIT (unitary sparse variant)
            mu_fine     = covguided.unitaryEspritSparse(Yb_fine, SoICols_fine, d, M);
            mu_wideFine = covguided.unitaryEspritSparse(Yb_wide, SoICols_wide, d, M);

            % Covariances for fine power estimation
            RY_fine     = (Yb_fine * Yb_fine') / N;     RY_fine     = (RY_fine + RY_fine')/2;
            RY_wideFine = (Yb_wide * Yb_wide') / N;     RY_wideFine = (RY_wideFine + RY_wideFine')/2;

            % Fine power estimation
            [~, ~, mu_fineSort,     ~] = covguided.estimatePowersFine(RY_fine,     SoICols_fine, mu_fine,     M);
            [~, ~, mu_wideFineSort, ~] = covguided.estimatePowersFine(RY_wideFine, SoICols_wide, mu_wideFine, M);

            % Per-source mu-errors (wrapped to [-pi,pi))
            err_mu_fine     = angle(exp(1i*(mu_true.' - mu_fineSort(:).')));      % 1xd
            err_mu_sector   = angle(exp(1i*(mu_true.' - mu_wideFineSort(:).')));  % 1xd
            
            %
            trial_fail_fine   = any(abs(err_mu_fine)   > thr_src);
            trial_fail_sector = any(abs(err_mu_sector) > thr_src);
            
            % Optional: keep RMS for MSE plots (but failures use per-source rule)
            e_vec = [ sqrt(mean(err_mu_fine.^2)), sqrt(mean(err_mu_sector.^2)) ];
            
            % Single consistent writes
            e_trials(iter, :)       = e_vec;
            fail_flags_fine(iter)   = trial_fail_fine;
            fail_flags_sector(iter) = trial_fail_sector;

        end
        %% === Per-δ aggregation ============================================
        % Per-method MSE (since e_trials are RMS per iter)
        MSE(sIdx, :, t) = mean(e_trials.^2, 1).';

        % Failure rates = fraction of trials where any source failed
        Fail_vs_delta(sIdx, 1, t) = mean(fail_flags_fine);
        Fail_vs_delta(sIdx, 2, t) = mean(fail_flags_sector);

        %% === Progress report =============================================
        ticksDone = ticksDone + 1;
        elapsed   = toc(tStart);
        eta       = elapsed / ticksDone * (totalTicks - ticksDone);
        fprintf('  [%2d/%2d] SNR=%+3d dB  delta=%+.4f rad  FR_cov=%5.1f%%  FR_sect=%5.1f%%  elapsed=%4ds  ETA=%4ds\n', ...
                ticksDone, totalTicks, ASNR_set(sIdx), delta(t), ...
                100*Fail_vs_delta(sIdx,1,t), 100*Fail_vs_delta(sIdx,2,t), ...
                round(elapsed), round(eta));
    end
end

%% === Global plotting style (IEEE-ish, single-column) =====================
fontName       = 'Times New Roman';  % or same as your main figs
fontSizeAxes   = 8;                  % tick labels
fontSizeLabel  = 9;                  % axis labels
fontSizeLegend = 8;
lineWidth      = 1.0;
axisLineWidth  = 0.75;
markerSize     = 4;

singleColWidth = 8.8;    % cm, IEEE single column ~8.8 cm
figHeight      = 6.5;    % cm, adjust as you like

% Color order (if you want explicit colors; otherwise MATLAB default is fine)
colors = lines(3);       % 1: cov-guided, 2: sector, 3: CRB


%% ========= RMSE vs delta (main A2 metric) =============
fig_rmse = figure('Color','w', ...
                  'Units','centimeters', ...
                  'Visible','off');
fig_rmse.Position(3:4) = [singleColWidth, figHeight];
set(fig_rmse, 'Renderer','painters');

tlo = tiledlayout(1, numel(ASNR_set), ...
                  'TileSpacing','compact', ...
                  'Padding','compact');

x = delta / (w/2);   % normalized offset

% legend handles + first axes
h_cov = []; h_sect = []; h_crb = [];
ax_first = [];

for sIdx = 1:numel(ASNR_set)
    ax = nexttile; 
    if sIdx == 1
        ax_first = ax;                     % remember first axes
    end
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    
    rmse_cov  = sqrt(squeeze(MSE(sIdx,1,:)));   % cov-guided
    rmse_sect = sqrt(squeeze(MSE(sIdx,2,:)));   % sectorization
    
    h1 = plot(ax, x, rmse_cov,  '-o', ...
        'LineWidth', lineWidth, ...
        'MarkerSize', markerSize, ...
        'Color', colors(1,:), ...
        'DisplayName','Cov-guided');
    
    h2 = plot(ax, x, rmse_sect, '-s', ...
        'LineWidth', lineWidth, ...
        'MarkerSize', markerSize, ...
        'Color', colors(2,:), ...
        'DisplayName','Sectorization');
    
    h3 = plot(ax, x, sqrt_CRB(sIdx,:), '-x', ...
        'LineWidth', lineWidth, ...
        'MarkerSize', markerSize, ...
        'Color', colors(3,:), ...
        'DisplayName','CRB');

    if sIdx == 1
        % remember handles for legend
        h_cov  = h1;
        h_sect = h2;
        h_crb  = h3;
    end

    xlabel(ax, '$\delta / (w_{\mathrm{sec}}/2)$', ...
           'FontName', fontName, ...
           'FontSize', fontSizeLabel, ...
           'Interpreter','latex');
    ylabel(ax, 'RMSE [rad]', ...
           'FontName', fontName, ...
           'FontSize', fontSizeLabel, ...
           'Interpreter','latex');

    title(ax, sprintf('ASNR = %d dB', ASNR_set(sIdx)), ...
          'FontName', fontName, ...
          'FontSize', fontSizeLabel, ...
          'Interpreter','latex');

    set(ax, 'FontName', fontName, ...
            'FontSize', fontSizeAxes, ...
            'LineWidth', axisLineWidth, ...
            'TickLabelInterpreter','latex');
end

% === Single compact legend for RMSE figure (first subplot only) =========
lgd = legend(ax_first, [h_cov, h_sect, h_crb], ...
    {'Prop Cov-guide', '[1] Sector', '$\sqrt{\mathrm{CRB}}$'}, ...
    'Location',    'southeast', ...          % inside top-right of first subplot
    'Orientation', 'vertical', ...
    'Box',         'on', ...
    'FontName',    fontName, ...
    'FontSize',    fontSizeLegend-1, ...     % a bit smaller than axes font
    'Interpreter', 'latex');

% Make tokens smaller so the legend is narrow
if isprop(lgd, 'ItemTokenSize')
    lgd.ItemTokenSize = [10 8];             % much smaller than default
end

drawnow;
exportgraphics(fig_rmse, fullfile(figDir,'A2_RMSE_vs_delta.pdf'), 'ContentType','vector');
close(fig_rmse);


%% ============ Failure rate vs delta ==============
fig_fail = figure('Color','w', ...
                  'Units','centimeters', ...
                  'Visible','off');
fig_fail.Position(3:4) = [singleColWidth, figHeight];
set(fig_fail, 'Renderer','painters');

tlo2 = tiledlayout(1, numel(ASNR_set), ...
                   'TileSpacing','compact', ...
                   'Padding','compact');

x = delta / (w/2);

h_cov = []; h_sect = [];
ax_first = [];

for sIdx = 1:numel(ASNR_set)
    ax = nexttile;
    if sIdx == 1
        ax_first = ax;
    end
    hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');
    
    fr_cov  = 100 * squeeze(Fail_vs_delta(sIdx,1,:));
    fr_sect = 100 * squeeze(Fail_vs_delta(sIdx,2,:));
    
    h1 = plot(ax, x, fr_cov,  '-o', ...
        'LineWidth', lineWidth, ...
        'MarkerSize', markerSize, ...
        'Color', colors(1,:), ...
        'DisplayName','Cov-guided');
    
    h2 = plot(ax, x, fr_sect, '-s', ...
        'LineWidth', lineWidth, ...
        'MarkerSize', markerSize, ...
        'Color', colors(2,:), ...
        'DisplayName','Sectorization');
    
    if sIdx == 1
        h_cov  = h1;
        h_sect = h2;
    end

    xlabel(ax, '$\delta / (w_{\mathrm{sec}}/2)$', ...
           'FontName', fontName, ...
           'FontSize', fontSizeLabel, ...
           'Interpreter','latex');
    ylabel(ax, 'Failure rate [\%]', ...
           'FontName', fontName, ...
           'FontSize', fontSizeLabel, ...
           'Interpreter','latex');
    
    title(ax, sprintf('ASNR = %d dB', ASNR_set(sIdx)), ...
          'FontName', fontName, ...
          'FontSize', fontSizeLabel, ...
          'Interpreter','latex');
    
    set(ax, 'FontName', fontName, ...
            'FontSize', fontSizeAxes, ...
            'LineWidth', axisLineWidth, ...
            'TickLabelInterpreter','latex');
end

% === Single compact legend for Failure-rate figure (first subplot) ======
lgd2 = legend(ax_first, [h_cov, h_sect], ...
    {'Prop Cov-guide', '[1] Sector'}, ...
    'Location',    'northeast', ...
    'Orientation', 'vertical', ...
    'Box',         'on', ...
    'FontName',    fontName, ...
    'FontSize',    fontSizeLegend-1, ...
    'Interpreter', 'latex');

if isprop(lgd2, 'ItemTokenSize')
    lgd2.ItemTokenSize = [10 8];
end

drawnow;
exportgraphics(fig_fail, fullfile(figDir,'A2_FailRate_vs_delta.pdf'), 'ContentType','vector');
close(fig_fail);


%% === Save results to CSV (A2) ============================================
% Long-format table: one row per (SNR, method, delta)
methodNames = ["Prop CovGuided","[11] Sector"];
delta_norm  = delta / (w/2);  % normalized x-axis used in plots

nrows = numel(ASNR_set) * numMethods * K;
T = table( ...
    'Size', [nrows 17], ...
    'VariableTypes', {'string','double','double','double','double','double','double', ...
                      'double','double','double','double','double','double','double','double','double','double'}, ...
    'VariableNames', {'Method','SNR_dB','Delta_rad','Delta_norm','FailureRate_pct','MSE','ThresholdRMS', ...
                      'M','N','NRF_coarse','NRF_fine','W','K','TotalIter','d','BoundaryIndex','SectorWidth_beams'});

r = 0;
for sIdx = 1:numel(ASNR_set)
    for m = 1:numMethods
        for t = 1:K
            r = r + 1;
            T.Method(r)            = methodNames(m);
            T.SNR_dB(r)            = ASNR_set(sIdx);
            T.Delta_rad(r)         = delta(t);
            T.Delta_norm(r)        = delta_norm(t);
            T.FailureRate_pct(r)   = 100 * Fail_vs_delta(sIdx, m, t);
            T.MSE(r)               = MSE(sIdx, m, t);
            T.ThresholdRMS(r)      = thrVal(sIdx, t);

            % config (replicated so the CSV is self-contained)
            T.M(r)                 = M;
            T.N(r)                 = N;
            T.NRF_coarse(r)        = NRF_coarse;
            T.NRF_fine(r)          = NRF_fine;
            T.W(r)                 = W;
            T.K(r)                 = K;
            T.TotalIter(r)         = totalIter;
            T.d(r)                 = d;
            T.BoundaryIndex(r)     = k_bdy;
            T.SectorWidth_beams(r) = W;
        end
    end
end

% Timestamped filename to avoid accidental overwrite
ts     = datestr(now, 'yyyymmdd_HHMMSS');
outCsv = fullfile(csvDir, sprintf('A2_sector_edge_results_M%d_W%d_%s.csv', M, W, ts));
writetable(T, outCsv);

fprintf('A2 results written to %s\n', outCsv);



%% ========================= HELPERS =========================
function idxOut = pad_truncate_to_nrf(idxIn, NRF_target, M)
% Ensure a contiguous block of NRF_target beams around the core set idxIn.
% If |idxIn|<NRF_target, pad contiguously; if >, take the tightest block.
    idxIn = unique(idxIn(:).', 'stable');
    if isempty(idxIn)
        idxOut = 1:NRF_target;
        return;
    end
    if numel(idxIn) >= NRF_target
        % take a contiguous block centered on the middle of idxIn
        mid  = idxIn(round(numel(idxIn)/2));
        half = floor(NRF_target/2);
        startIdx = mid - half;
        idxOut = mod((startIdx:(startIdx+NRF_target-1)) - 1, M) + 1;
        return;
    end
    % padding needed
    need = NRF_target - numel(idxIn);
    lo = min(idxIn); hi = max(idxIn);
    padCandidates = setdiff(mod((lo-need):(hi+need) - 1, M) + 1, idxIn, 'stable');
    idxOut = [idxIn, padCandidates(1:need)];
end
