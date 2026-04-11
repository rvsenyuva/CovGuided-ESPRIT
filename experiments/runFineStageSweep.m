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

%% ------- Invariants (compute once) -------
M   = 32;
N   = 100;
NRF = 12;

%% --- Fine-stage beam budget sweep (S1) ---
% candidate fine-stage beam budgets
kf_list = [2 3 4];
numKf   = numel(kf_list);

%%
opts_sel = struct('eig_tol_scale',1e-6, 'tikhonov_gamma', 0, ...
                  'alpha',0.5, 'normalize_score', true, ...
                  'prune_topq', 0);                      % 0 = off; try 3 for ~3–5× fewer windows

%% RF mask (element-space subsampling for "coarse")
maskIdx  = covguided.centeredContiguousMask(M, NRF);
mIdxFlip = NRF:-1:1;                                     % for FBA reflect

%% Precompute once (outside parfor)
phase_k  = exp(1i*(M-1)*pi*(0:M-1).'/M);                 % size Mx1, for row-wise scaling

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

%% "Oracle" beam set (invariant)
beamGroups_ini  = covguided.sectorizeDftBeams(mu_true, M, 4, 1);
SoICols_ini     = covguided.selectAdjacentPairsFromSectorizedCov(R_signal, beamGroups_ini);
SoICols_true    = unique(cell2mat(SoICols_ini)).';

%% SNR & noise
ASNR     = [0 5 10 15];
noiseVar = real(trace(R_signal)) ./ (M * 10.^(ASNR/10));

%% CRB
CRB = zeros(1,length(noiseVar));
PorthA = eye(M)-Atrue*pinv(Atrue);
derA = (1i*transpose(0:M-1)).*Atrue;
H = transpose(derA'*PorthA*derA);
for noiseInd=1:length(noiseVar)
    secondTerm = ((Atrue'*Atrue)*R_source)/(noiseVar(noiseInd));
    firstTerm  = eye(d)+secondTerm;
    CRBinvterm = diag(inv(real((R_source*(firstTerm\secondTerm)).*H)));
    CRB(noiseInd) = sum((noiseVar(noiseInd)/(2*N))*CRBinvterm)/d;
end
sqrtCRB = sqrt(CRB);

%% --- Accumulators -------------------------------------------------------
totalIter     = 1e4;

MSE_fine      = zeros(numKf, length(ASNR));       % RMSE^2 per k_f and SNR
MSE_wide      = zeros(numKf, length(ASNR));       % RMSE^2 per k_f and SNR
failRate_fine = nan(numKf, length(ASNR));         % failure probability per k_f and SNR
failRate_wide = nan(numKf, length(ASNR));         % failure probability per k_f and SNR
ERR_fine      = cell(numKf, length(ASNR));        % store per-trial RMS errors if desired
ERR_wide      = cell(numKf, length(ASNR));        % store per-trial RMS errors if desired

kThr          = 3;

%% --- Outer SNR loop -----------------------------------------------------
for snrIndex = 1:length(ASNR)

    %% --- S1: per-trial containers for this SNR ---
    e_trials_fine  = nan(totalIter, numKf);           % per trial, per k_f
    e_trials_wide  = nan(totalIter, numKf);
    nvStd          = sqrt(noiseVar(snrIndex)/2); 

    %% --- Setup progress bar for this SNR run ---
    nIter = totalIter;
    currentASNR = ASNR(snrIndex);
    q = parallel.pool.DataQueue;
    
    afterEach(q, @(~) updateProgress(1, nIter, currentASNR));  % pass both nIter & ASNR
    updateProgress(NaN, nIter, currentASNR);                   % reset bar for this SNR


    %% PARFOR over Monte-Carlo iterations
    parfor iter = 1:totalIter
        % Make a stream whose *substream* encodes the iteration
        s = RandStream('Threefry','Seed', 5489);
        s.Substream = iter;                                                % 1,2,3,... ensures independence per iteration
        RandStream.setGlobalStream(s);                                     % activate for rand/randn/randi

        %% --- Coarse stage synthetic data (element space) ---
        S = Psqrt * (randn(d,N) + 1i*randn(d,N));
        X = Atrue * S;

        %% noise and received
        Z = nvStd * (randn(M,N) + 1i*randn(M,N));
        Y_coarse = X + Z;

        %% element-space subarray "beamforming"
        Yb_coarse = Y_coarse(maskIdx, :);           % NRF x N
        Yb_wideES = Y_coarse(1:NRF, :);             % baseline wide (as in your code)

        %% empirical covariances (Hermitian once)
        RY_coarse = (Yb_coarse * Yb_coarse') / N;   RY_coarse = (RY_coarse + RY_coarse') / 2;
        RY_wideES = (Yb_wideES * Yb_wideES') / N;   RY_wideES = (RY_wideES + RY_wideES') / 2;

        %% FBA (coarse)
        RY_fba    = (RY_coarse + conj(RY_coarse(mIdxFlip, mIdxFlip))) * 0.5;  % already Hermitian

        %% TLS-ESPRIT on coarse/wideES covariances
        mu_coarse = covguided.tlsEsprit(RY_fba,   d);
        mu_wide   = covguided.tlsEsprit(RY_wideES, d);

        %% Coarse power estimates (as in your code)
        [~, ~, mu_courseSort, R_course] = covguided.estimatePowersCoarse(RY_fba,    maskIdx,  mu_coarse, M);
        [~, ~, mu_wideSort,   R_wide]   = covguided.estimatePowersCoarse(RY_wideES, 1:NRF,     mu_wide,   M);

        %% --- Fine stage synthetic data (element space) ---
        Z2     = nvStd * (randn(M,N) + 1i*randn(M,N));
        S2     = Psqrt * (randn(d,N) + 1i*randn(d,N));
        X2     = Atrue * S2;
        Y_fine = X2 + Z2;

        %% --- FFT path for W_SE' * Y_fine (big win) ---
        Yb_full  = phase_k .* fft(Y_fine, M, 1);

        %%
        e_vec_fine = nan(1, numKf);
        e_vec_wide = nan(1, numKf);

        %% =========================================================
        %  S1: loop over fine-stage beam budgets k_f
        %  =========================================================
        for kk = 1:numKf
             %%
             kf = kf_list(kk);
             %% Contiguous k_f-beam sectors around each coarse DoA
             beamGroups_fine = covguided.sectorizeDftBeams(mu_courseSort, M, 2*kf, 1);
             selBeams_fine   = covguided.selectBeamsetsFromSectorizedCov(R_course, beamGroups_fine, kf, opts_sel);
             SoICols_fine    = unique(cell2mat(selBeams_fine)).';

             %% Contiguous k_f-beam sectors around each coarse DoA
             beamGroups_wide = covguided.sectorizeDftBeams(mu_wideSort, M, kf, 1);
             SoICols_wide    = unique(cell2mat(beamGroups_wide)).';

             %% Beamspace restriction for this k_f
             Yb_fine = Yb_full(SoICols_fine, :);
             Yb_wide = Yb_full(SoICols_wide, :);
      
             %% Unitary ESPRIT on fine and wide
             mu_fine = covguided.unitaryEspritSparse(Yb_fine, SoICols_fine, d, M);
             mu_wide = covguided.unitaryEspritSparse(Yb_wide, SoICols_wide, d, M);
            
             %% Covariances for power estimation (fine)
             RY_fine = (Yb_fine * Yb_fine') / N;      
             RY_fine = (RY_fine + RY_fine') / 2;

             %% Covariances for power estimation (wide)
             RY_wide = (Yb_wide * Yb_wide') / N;
             RY_wide = (RY_wide + RY_wide') / 2;

             %% Power estimation
             [~, ~, mu_fineSort, ~]      = covguided.estimatePowersFine(RY_fine, SoICols_fine, mu_fine, M);
             [~, ~, mu_wideSort_fine, ~] = covguided.estimatePowersFine(RY_wide, SoICols_wide, mu_wide, M);

             %% --- Per-trial RMS error and failure flag for this k_f ---
             e_vec_fine(kk) = sqrt(mean(angle(exp(1i*(mu_true - mu_fineSort))).^2));
             e_vec_wide(kk) = sqrt(mean(angle(exp(1i*(mu_true - mu_wideSort_fine))).^2));

        end
        % Store S1 results for this trial
        e_trials_fine(iter, :) = e_vec_fine;
        e_trials_wide(iter, :) = e_vec_wide;

        % update progress bar
        send(q, 1);   
    end % parfor

    %% ==============================================================
    %  Reduce over trials: RMSE^2 and failure rate vs k_f
    % ==============================================================
    MSE_fine(:, snrIndex)      = mean(e_trials_fine.^2, 1).';
    MSE_wide(:, snrIndex)      = mean(e_trials_wide.^2, 1).';

    thr                        = kThr * max(sqrtCRB(snrIndex), 1e-12);

    failRate_fine(:, snrIndex) = mean(e_trials_fine > thr, 1).';
    failRate_wide(:, snrIndex) = mean(e_trials_wide > thr, 1).';

    % If you also want per-trial error vectors stored:
    for kk = 1:numKf
        ERR_fine{kk, snrIndex} = e_trials_fine(:, kk);
        ERR_wide{kk, snrIndex} = e_trials_wide(:, kk);
    end
    
end % SNR loop

%% ================== S1 PLOTTING (place BEFORE updateProgress) ==================
% Global plotting style (consistent with A2 script)
fontName       = 'Times New Roman';
fontSizeAxes   = 8;          % tick labels
fontSizeLabel  = 9;          % axis labels
fontSizeLegend = 8;
lineWidth      = 1.0;
axisLineWidth  = 0.75;
markerSize     = 4;

legendTokenSize_rmse  = [4 2];    % was [10 8], shorter legend symbols
legendTokenSize_fail  = [6 4];    % was [10 8], shorter legend symbols

singleColWidth = 8.8;        % cm, IEEE single column
figHeight      = 6.5;        % cm

% Colors: one per K_f; CRB uses its own color
colors = lines(numKf);       % color(k,:) used for K_f = kf_list(k)

RMSE_fine = sqrt(MSE_fine);  % [numKf x length(ASNR)]
RMSE_wide = sqrt(MSE_wide);  % [numKf x length(ASNR)]


%% ========= S1: RMSE vs ASNR (single-column figure, both methods) =======
fig_rmse_S1 = figure('Color','w', ...
                     'Units','centimeters', ...
                     'Visible','off');
fig_rmse_S1.Position(3:4) = [singleColWidth, figHeight];
set(fig_rmse_S1, 'Renderer','painters');

ax = axes('Parent', fig_rmse_S1);
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

x = ASNR;   % ASNR in dB

h_fine = gobjects(numKf,1);
h_wide = gobjects(numKf,1);

for kk = 1:numKf
    % Proposed covariance-guided selector
    h_fine(kk) = plot(ax, x, RMSE_fine(kk,:), '-o', ...
        'LineWidth',  lineWidth, ...
        'MarkerSize', markerSize, ...
        'Color',      colors(kk,:), ...
        'DisplayName', sprintf('Cov-guided, $K_f = %d$', kf_list(kk)));

    % Pure sectorization baseline
    h_wide(kk) = plot(ax, x, RMSE_wide(kk,:), '--s', ...
        'LineWidth',  lineWidth, ...
        'MarkerSize', markerSize, ...
        'Color',      colors(kk,:), ...
        'DisplayName', sprintf('Sector, $K_f = %d$', kf_list(kk)));
end

% CRB curve (grey / black)
h_crb = plot(ax, x, sqrtCRB, '-x', ...
    'LineWidth',  lineWidth, ...
    'MarkerSize', markerSize, ...
    'Color',      [0 0 0], ...
    'DisplayName','CRB');

xlabel(ax, 'ASNR [dB]', ...
       'FontName', fontName, ...
       'FontSize', fontSizeLabel, ...
       'Interpreter','latex');

ylabel(ax, 'RMSE [rad]', ...
       'FontName', fontName, ...
       'FontSize', fontSizeLabel, ...
       'Interpreter','latex');

%title(ax, 'S1: RMSE vs ASNR', ...
%      'FontName', fontName, ...
%      'FontSize', fontSizeLabel, ...
%      'Interpreter','latex');

set(ax, 'FontName', fontName, ...
        'FontSize', fontSizeAxes, ...
        'LineWidth', axisLineWidth, ...
        'TickLabelInterpreter','latex');

set(ax, 'YScale','log');        % log scale on RMSE
% Tweak limits as needed after you see the data:
ylim(ax, [1e-3 2]);             
yticks(ax, [1e-3 1e-2 1e-1 1]);

lgd_S1_rmse = legend(ax, [h_fine(:); h_wide(:); h_crb], ...
    [arrayfun(@(k)sprintf('Prop Cov-guide, $K_f = %d$', k), kf_list, 'UniformOutput',false), ...
     arrayfun(@(k)sprintf('[1] Sector, $K_f = %d$', k), kf_list, 'UniformOutput',false), ...
     {'$\sqrt{\mathrm{CRB}}$'}], ...
    'Location',    'northeast', ...
    'Orientation', 'vertical', ...
    'Box',         'on', ...
    'FontName',    fontName, ...
    'FontSize',    fontSizeLegend, ...
    'Interpreter', 'latex');

if isprop(lgd_S1_rmse, 'ItemTokenSize')
    lgd_S1_rmse.ItemTokenSize = legendTokenSize_rmse;
end

drawnow;
exportgraphics(fig_rmse_S1, fullfile(figDir,'S1_RMSE_vs_ASNR.pdf'), 'ContentType','vector');
close(fig_rmse_S1);

%% ============ S1: Failure rate vs ASNR (single-column, both methods) ====
fig_fail_S1 = figure('Color','w', ...
                     'Units','centimeters', ...
                     'Visible','off');
fig_fail_S1.Position(3:4) = [singleColWidth, figHeight];
set(fig_fail_S1, 'Renderer','painters');

ax2 = axes('Parent', fig_fail_S1);
hold(ax2, 'on'); grid(ax2, 'on'); box(ax2, 'on');

x = ASNR;   % ASNR in dB

h_fine_fail = gobjects(numKf,1);
h_wide_fail = gobjects(numKf,1);

for kk = 1:numKf
    fr_fine = 100 * failRate_fine(kk,:);   % in percent
    fr_wide = 100 * failRate_wide(kk,:);   % in percent

    h_fine_fail(kk) = plot(ax2, x, fr_fine, '-o', ...
        'LineWidth',  lineWidth, ...
        'MarkerSize', markerSize, ...
        'Color',      colors(kk,:), ...
        'DisplayName', sprintf('Cov-guided, $K_f = %d$', kf_list(kk)));

    h_wide_fail(kk) = plot(ax2, x, fr_wide, '--s', ...
        'LineWidth',  lineWidth, ...
        'MarkerSize', markerSize, ...
        'Color',      colors(kk,:), ...
        'DisplayName', sprintf('Sector, $K_f = %d$', kf_list(kk)));
end

xlabel(ax2, 'ASNR [dB]', ...
       'FontName', fontName, ...
       'FontSize', fontSizeLabel, ...
       'Interpreter','latex');

ylabel(ax2, 'Failure rate [\%]', ...
       'FontName', fontName, ...
       'FontSize', fontSizeLabel, ...
       'Interpreter','latex');

%title(ax2, 'S1: Failure rate vs ASNR', ...
%      'FontName', fontName, ...
%      'FontSize', fontSizeLabel, ...
%      'Interpreter','latex');

set(ax2, 'FontName', fontName, ...
         'FontSize', fontSizeAxes, ...
         'LineWidth', axisLineWidth, ...
         'TickLabelInterpreter','latex');

% You can switch to log scale if you prefer:
% set(ax2, 'YScale', 'log');

lgd_S1_fail = legend(ax2, [h_fine_fail(:); h_wide_fail(:)], ...
    [arrayfun(@(k)sprintf('Prop Cov-guide, $K_f = %d$', k), kf_list, 'UniformOutput',false), ...
     arrayfun(@(k)sprintf('[1] Sector, $K_f = %d$', k), kf_list, 'UniformOutput',false)], ...
    'Location',    'northeast', ...
    'Orientation', 'vertical', ...
    'Box',         'on', ...
    'FontName',    fontName, ...
    'FontSize',    fontSizeLegend, ...
    'Interpreter', 'latex');

if isprop(lgd_S1_fail, 'ItemTokenSize')
    lgd_S1_fail.ItemTokenSize = legendTokenSize_fail;  % [6 4]
end

drawnow;
exportgraphics(fig_fail_S1, fullfile(figDir,'S1_FailRate_vs_ASNR.pdf'), 'ContentType','vector');
close(fig_fail_S1);

%% ================== SAVE S1 RESULTS TO CSV ============================
% ASNR        : 1 x S
% kf_list     : 1 x numKf
% MSE_fine    : numKf x S
% MSE_wide    : numKf x S
% failRate_*  : numKf x S
% sqrtCRB     : 1 x S
% totalIter   : scalar
% kThr        : scalar
%
% Creates:
%   - S1_aggregates.csv   (one row per {method, Kf, ASNR})
%   - S1_trials_subset.csv (per-trial errors for selected ASNRs, both methods)

% ---------- Aggregated metrics per {method, Kf, ASNR} ----------
[numKf, numASNR] = size(MSE_fine);
[KK, AA] = ndgrid(kf_list(:), ASNR(:));          % both [numKf x numASNR]

RMSE_fine = sqrt(MSE_fine);                      % [numKf x numASNR]
RMSE_wide = sqrt(MSE_wide);                      % [numKf x numASNR]
sqrtCRB_rep = repmat(sqrtCRB(:).', numKf, 1);    % same size as MSE_*

% Cov-guided entries
Agg = table;
Agg.Kf         = KK(:);
Agg.ASNR_dB    = AA(:);
Agg.method     = repmat("Cov-guided", numKf*numASNR, 1);
Agg.RMSE_rad   = RMSE_fine(:);
Agg.FailRate   = failRate_fine(:);
Agg.sqrtCRB_rd = sqrtCRB_rep(:);
Agg.kThr       = kThr * ones(numKf*numASNR,1);
Agg.totalIter  = totalIter * ones(numKf*numASNR,1);

% Sectorization entries
Agg2 = table;
Agg2.Kf         = KK(:);
Agg2.ASNR_dB    = AA(:);
Agg2.method     = repmat("Sectorization", numKf*numASNR, 1);
Agg2.RMSE_rad   = RMSE_wide(:);
Agg2.FailRate   = failRate_wide(:);
Agg2.sqrtCRB_rd = sqrtCRB_rep(:);
Agg2.kThr       = kThr * ones(numKf*numASNR,1);
Agg2.totalIter  = totalIter * ones(numKf*numASNR,1);

Agg = [Agg; Agg2];

writetable(Agg, fullfile(csvDir,'S1_aggregates.csv'));

% ---------- Optional: per-trial subset for ECDF / scatter -------------
wantASNR = [0 15];   % adjust as you like

Trials = table;
for a = wantASNR
    s = find(ASNR == a, 1);
    if isempty(s)
        warning('ASNR %g dB not found in ASNR vector, skipping.', a);
        continue;
    end
    for kk = 1:numKf
        e_f = ERR_fine{kk, s};   % cov-guided
        e_w = ERR_wide{kk, s};   % sectorization

        if ~isempty(e_f)
            ntr = numel(e_f);
            T_f = table( ...
                repmat(a, ntr, 1), ...
                repmat(kf_list(kk), ntr, 1), ...
                (1:ntr).', ...
                e_f(:), ...
                repmat("Cov-guided", ntr, 1), ...
                'VariableNames', {'ASNR_dB','Kf','trial','err_rms_rad','method'});
            Trials = [Trials; T_f]; %#ok<AGROW>
        end

        if ~isempty(e_w)
            ntr = numel(e_w);
            T_w = table( ...
                repmat(a, ntr, 1), ...
                repmat(kf_list(kk), ntr, 1), ...
                (1:ntr).', ...
                e_w(:), ...
                repmat("Sectorization", ntr, 1), ...
                'VariableNames', {'ASNR_dB','Kf','trial','err_rms_rad','method'});
            Trials = [Trials; T_w]; %#ok<AGROW>
        end
    end
end

if ~isempty(Trials)
    writetable(Trials, fullfile(csvDir,'S1_trials_subset.csv'));
end

%% ========================================================================
% Local progress bar function with ASNR label (works for scripts)
% ========================================================================
function updateProgress(incr, nIter, asnrVal)

persistent count lastPct
if isnan(incr)
    count = 0; lastPct = 0;
    barWidth = 40;
    fprintf('\rASNR = %3.0f dB  [%s]   0%%', asnrVal, repmat(' ',1,barWidth));
    return
end

if isempty(count), count = 0; lastPct = 0; end
count = count + 1;

pct = floor(100 * count / nIter);
if pct >= lastPct + 5 || count == nIter
    lastPct = pct;
    barWidth = 40;
    nBars = round(pct/100 * barWidth);
    fprintf('\rASNR = %3.0f dB  [%s%s] %3d%%', asnrVal, ...
        repmat('=',1,nBars), repmat(' ',1,barWidth-nBars), pct);
    if count == nIter
        fprintf('\n');  % move to new line at completion
    end
end

end