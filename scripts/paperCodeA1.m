clearvars;
%% --- Invariants (compute once) ------------------------------------------
M   = 32;
N   = 100;
NRF = 12;

%% set to 2 for 12->6; increase to {3,4,...} for ablations (S1)
K_FINE   = [2 4];
%opts_sel = struct('eig_tol_scale', 1e-6, 'tikhonov_gamma', 0, 'alpha', 0.5, 'normalize_score', true);
opts_sel = struct('eig_tol_scale',1e-6, 'tikhonov_gamma',0, ...
                  'alpha',0.5, 'normalize_score', true, ...
                  'prune_topq', 4);          % 0 = off; try 3 for ~3–5× fewer windows

%% RF mask (element-space subsampling for "coarse")
maskIdx  = centered_contiguous_mask(M, NRF);
mIdxFlip = NRF:-1:1;                                      % for FBA reflect

%% Precompute once (outside parfor)
phase_k  = exp(1i*(M-1)*pi*(0:M-1).'/M);                  % size Mx1, for row-wise scaling

%% Source powers, true DOAs
pows     = [0.95, 0.5, 0.1];     
R_source = diag(pows);
Psqrt    = diag(sqrt(pows/2));
mu_true  = [-2.1; 0.5; 2.5];
d        = numel(mu_true);

%% Steering (precompute once)
Afun     = @(mu) exp(1i * (0:M-1).' * mu.'); 
Atrue    = Afun(mu_true);

%% Signal covariance (invariant)
R_signal = Atrue * R_source * Atrue';

%% SNR & noise
ASNR     = -4:1:6;
noiseVar = real(trace(R_signal))./(M*10.^(ASNR/10));

%% --- Accumulators -------------------------------------------------------
totalIter   = 1e4;
numMethods  = 4;

MSE         = zeros(numMethods, length(ASNR), length(K_FINE));

med_t_cov   = zeros(2, length(ASNR), length(K_FINE));
med_t_sel   = zeros(2, length(ASNR), length(K_FINE));
med_t_es    = zeros(2, length(ASNR), length(K_FINE));
med_t_total = zeros(2, length(ASNR), length(K_FINE)); 

%% --- Outer SNR loop -----------------------------------------------------
for k = 1:length(K_FINE)

    K_FINE_CURRENT = K_FINE(k);

    for snrIndex = 1:length(ASNR)

        % Per-SNR local accumulators (parfor-friendly)
        e_trials    = zeros(totalIter, numMethods);
        nvStd       = sqrt(noiseVar(snrIndex)/2); 

        % Timing arrays (per-trial inside parfor)
        t_cov_trials   = zeros(totalIter, 2);
        t_sel_trials   = zeros(totalIter, 2);
        t_es_trials    = zeros(totalIter, 2);
        t_total_trials = zeros(totalIter, 2);

        % --- Setup progress bar for this SNR run ---
        nIter = totalIter;
        currentASNR = ASNR(snrIndex);
        q = parallel.pool.DataQueue;
    
        afterEach(q, @(~) updateProgress(1, nIter, currentASNR));  % pass both nIter & ASNR
        updateProgress(NaN, nIter, currentASNR);                   % reset bar for this SNR

        %%% PARFOR over Monte-Carlo iterations
        parfor iter = 1:totalIter
            % Make a stream whose *substream* encodes the iteration
            s = RandStream('Threefry','Seed', 5489);
            s.Substream = iter;                 % 1,2,3,... ensures independence per iteration
            RandStream.setGlobalStream(s);      % activate for rand/randn/randi

            %%% --- Coarse stage synthetic data (element space) ---
            S = Psqrt * (randn(d,N) + 1i*randn(d,N));
            X = Atrue * S;

            % noise and received
            Z = nvStd * (randn(M,N) + 1i*randn(M,N));
            Y_coarse = X + Z;

            % element-space subarray "beamforming"
            Yb_coarse = Y_coarse(maskIdx, :);           % NRF x N
            Yb_wideES = Y_coarse(1:NRF, :);             % baseline wide (as in your code)

            % empirical covariances (Hermitian once)
            RY_coarse = (Yb_coarse * Yb_coarse') / N;   RY_coarse = (RY_coarse + RY_coarse')/2;
            RY_wideES = (Yb_wideES * Yb_wideES') / N;   RY_wideES = (RY_wideES + RY_wideES')/2;

            % FBA (coarse)
            RY_fba = (RY_coarse + conj(RY_coarse(mIdxFlip, mIdxFlip))) * 0.5;  % already Hermitian

            % TLS-ESPRIT on coarse/wideES covariances
            mu_coarse = tls_esprit(RY_fba,   d);
            mu_wide   = tls_esprit(RY_wideES, d);

            % Coarse power estimates
            t_cov_tic1                       = tic;
            [~, ~, mu_courseSort, R_course]  = estimate_powers_course(RY_fba,    maskIdx,  mu_coarse, M);
            t_cov1                           = toc(t_cov_tic1);
        
            % Coarse power estimates
            t_cov_tic2                     = tic;
            [~, ~, mu_wideSort,   ~]       = estimate_powers_course(RY_wideES, 1:NRF,     mu_wide,   M);
            t_cov2                         = toc(t_cov_tic2);
    
            % covariance-based beam selection for fine stage
            t_sel_tic       = tic;
            beamGroups_fine = sectorize_dft_beams(mu_courseSort, M, 4, 1);
            selBeams_fine   = select_beamsets_from_sectorized_cov(R_course, beamGroups_fine, K_FINE_CURRENT, opts_sel);
            SoICols_fine    = unique(cell2mat(selBeams_fine)).';
            t_sel           = toc(t_sel_tic);

            % sectorization for fine stage
            t_sec_tic       = tic;
            beamGroups_wide = sectorize_dft_beams(mu_wideSort,   M, K_FINE_CURRENT, 1);
            SoICols_wide    = unique(cell2mat(beamGroups_wide)).';
            t_sec           = toc(t_sec_tic);
    
            %%% --- Fine stage synthetic data (element space) ---
            Z2     = nvStd * (randn(M,N) + 1i*randn(M,N));
            S2     = Psqrt * (randn(d,N) + 1i*randn(d,N));
            X2     = Atrue * S2;
            Y_fine = X2 + Z2;
    
            %%% --- FFT path for W_SE' * Y_fine (big win) ---
            Yb_full = phase_k .* fft(Y_fine, M, 1);
    
            Yb_fine = Yb_full(SoICols_fine, :);
            Yb_wide = Yb_full(SoICols_wide, :);

            % Unitary ESPRIT on fine
            t_es_tic1    = tic;
            mu_fine     = unitary_esprit_sparse(Yb_fine, SoICols_fine, d, M);
            t_es1        = toc(t_es_tic1);
    
            % Unitary ESPRIT on wide beams
            t_es_tic2   = tic;
            mu_wideFine = unitary_esprit_sparse(Yb_wide, SoICols_wide, d, M);
            t_es2       = toc(t_es_tic2);
    
            % Covariances for power estimation (fine)
            RY_fine      = (Yb_fine * Yb_fine') / N;      RY_fine      = (RY_fine      + RY_fine')/2;
            RY_wideFine  = (Yb_wide * Yb_wide') / N;      RY_wideFine  = (RY_wideFine  + RY_wideFine')/2;
    
            % power estimation
            [~, ~, mu_fineSort,     R_fine]     = estimate_powers_fine(RY_fine,     SoICols_fine, mu_fine,     M);
            [~, ~, mu_wideFineSort, R_wideFine] = estimate_powers_fine(RY_wideFine, SoICols_wide, mu_wideFine, M);
            
            % --- build local row results (scalars -> 1x6 vectors)
            e_vec = [ ...
            sqrt(mean(angle(exp(1i*(mu_true - mu_courseSort   ))).^2)), ...
            sqrt(mean(angle(exp(1i*(mu_true - mu_wideSort     ))).^2)), ...
            sqrt(mean(angle(exp(1i*(mu_true - mu_fineSort     ))).^2)), ...
            sqrt(mean(angle(exp(1i*(mu_true - mu_wideFineSort ))).^2)), ...
            ];
    
            % --- single consistent sliced writes
            e_trials(iter,   :)    = e_vec;
    
            % bookkeeping
            t_cov_trials(iter,:)   = [t_cov1 t_cov2];
            t_sel_trials(iter,:)   = [t_sel t_sec];
            t_es_trials(iter,:)    = [t_es1 t_es2];
            t_total_trials(iter,:) = [t_cov1 + t_sel + t_es1 t_cov2 + t_sec + t_es2];
    
            % update progress bar
            send(q, 1);   

        end % parfor

        %%% --- Reduce after parfor (vectorized) ---
        % RMSE^2 (MSE) per method at this ASNR
        MSE(:, snrIndex, k)            = mean(e_trials.^2, 1).';  % since e_trials is already RMS per iter
    
        % --- Timing medians (ms) for Pareto A1 ---
        for m = 1:2
            med_t_cov(m, snrIndex, k)   = 1e3 * median(t_cov_trials(:,m));
            med_t_sel(m, snrIndex, k)   = 1e3 * median(t_sel_trials(:,m));
            med_t_es(m, snrIndex, k)    = 1e3 * median(t_es_trials(:,m));
            med_t_total(m, snrIndex, k) = 1e3 * median(t_total_trials(:,m));
        end
    
    end % SNR loop
end
%% --- Post-processing & plots (unchanged except cosmetic) ----------------
RMSE    = sqrt(MSE);

%% Pareto A1
% Choose the two SNRs to show (must exist in your ASNR vector)
snr_use = [3 6];
sIdx    = arrayfun(@(x) find(ASNR==x, 1, 'first'), snr_use);

% Assemble [C x S] matrices in the order of condLabels
condLabels = {'Cov 12->6','Cov 12->12','Sect 12->6','Sect 12->12'};

med_t_total_ms = [
    med_t_total(1, sIdx, 1);
    med_t_total(1, sIdx, 2);
    med_t_total(2, sIdx, 1);
    med_t_total(2, sIdx, 2);
   ];

rmse_rad = [
    RMSE(3, sIdx, 1);
    RMSE(3, sIdx, 2);
    RMSE(4, sIdx, 1);
    RMSE(4, sIdx, 2);
    ];

% Optional breakdown (same shape) for the side timing table
bk.cov_ms   = [
    med_t_cov(1, sIdx, 1);
    med_t_cov(1, sIdx, 2);
    med_t_cov(2, sIdx, 1);
    med_t_cov(2, sIdx, 2);
    ];
bk.sel_ms   = [
    med_t_sel(1, sIdx, 1);
    med_t_sel(1, sIdx, 2);
    med_t_sel(2, sIdx, 1);
    med_t_sel(2, sIdx, 2);
    ];
bk.es_ms    = [
    med_t_es(1, sIdx, 1);
    med_t_es(1, sIdx, 2);
    med_t_es(2, sIdx, 1);
    med_t_es(2, sIdx, 2);
    ];

bk.total_ms = med_t_total_ms;

% Sanity: print mapping from labels to rows
fprintf('--- Condition rows ---\n');
for c = 1:numel(condLabels)
    fprintf('%2d: %s\n', c, condLabels{c});
end

% Find row and column for "Sect 12->12 @ 3 dB"
rowSect1212 = find(strcmp(condLabels, 'Sect 12->12'));
col3dB      = find(snr_use == 3);

fprintf('\nSect 12->12 @ 3 dB:\n');
fprintf('  t_total_ms = %.4f\n', med_t_total_ms(rowSect1212, col3dB));
fprintf('  RMSE_rad   = %.6f\n', rmse_rad(rowSect1212, col3dB));


outdir = 'figures_pdf';
if ~exist(outdir,'dir'), mkdir(outdir); end

%% Style for A1 Pareto plot (match A2 figures)
plot_pareto_A1(med_t_total_ms, rmse_rad, ...
    'snr', snr_use, ...
    'condLabels', condLabels, ...
    'breakdown', bk, ...
    'connect_all_per_snr', true, ...  % dashed SNR-wise connectors
    'print_hull', true, ...
    'outpdf', fullfile('figures_pdf','A1_pareto.pdf'), ...
    'tablecsv', fullfile('figures_pdf','A1_timing_table.csv'));

% ========================================================================
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