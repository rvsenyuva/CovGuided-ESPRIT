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

%% --- Invariants (compute once) ------------------------------------------
M   = 32;
N   = 100;
NRF = 12;

%% RF mask (element-space subsampling for "coarse")
maskIdx = covguided.centeredContiguousMask(M, NRF);
mIdxFlip = NRF:-1:1;                                     % for FBA reflect

%% Precompute once (outside parfor)
phase_k = exp(1i*(M-1)*pi*(0:M-1).'/M);                  % size Mx1, for row-wise scaling

%% Source powers, true DOAs
pows     = [0.95, 0.5, 0.1];     
R_source = diag(pows);
Psqrt = diag(sqrt(pows/2));
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
ASNR     = -10:1:20;
noiseVar = real(trace(R_signal))./(M*10.^(ASNR/10));

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
totalIter  = 1e4;
numMethods = 5;

MSE        = zeros(numMethods, length(ASNR));
ERR        = cell(numMethods, length(ASNR));
LPA_TRIAL  = cell(numMethods, length(ASNR));
failRate   = nan(numMethods, length(ASNR)); 
k          = 3;

%% --- Outer SNR loop -----------------------------------------------------
for snrIndex = 1:length(ASNR)

    % Per-SNR local accumulators (parfor-friendly)
    e_trials    = zeros(totalIter, numMethods);
    lpa_trials  = zeros(totalIter, numMethods);
    nvStd       = sqrt(noiseVar(snrIndex)/2); 

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
        % signal
        S = Psqrt * (randn(d,N) + 1i*randn(d,N));
        X = Atrue * S;  % uses precomputed Atrue

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
        mu_coarse = covguided.tlsEsprit(RY_fba,   d);
        mu_wide   = covguided.tlsEsprit(RY_wideES, d);

        % Coarse power estimates (as in your code)
        [~, ~, mu_courseSort, R_course]     = covguided.estimatePowersCoarse(RY_fba,    maskIdx,  mu_coarse, M);
        [~, ~, mu_wideSort,   R_wide]       = covguided.estimatePowersCoarse(RY_wideES, 1:NRF,     mu_wide,   M);

        % Sectorization for fine stage
        beamGroups_fine = covguided.sectorizeDftBeams(mu_courseSort, M, 4, 1);
        selBeams_fine   = covguided.selectAdjacentPairsFromSectorizedCov(R_course, beamGroups_fine);
        SoICols_fine    = unique(cell2mat(selBeams_fine)).';

        beamGroups_wide = covguided.sectorizeDftBeams(mu_wideSort,   M, 2, 1);
        SoICols_wide    = unique(cell2mat(beamGroups_wide)).';

        %%% --- Fine stage synthetic data (element space) ---
        Z2 = nvStd * (randn(M,N) + 1i*randn(M,N));
        S2 = Psqrt * (randn(d,N) + 1i*randn(d,N));
        X2 = Atrue * S2;
        Y_fine = X2 + Z2;

        %%% --- FFT path for W_SE' * Y_fine (big win) ---
        Yb_full = phase_k .* fft(Y_fine, M, 1);

        Yb_fine = Yb_full(SoICols_fine, :);
        Yb_wide = Yb_full(SoICols_wide, :);
        Yb_true = Yb_full(SoICols_true, :);

        % Unitary ESPRIT on fine/wide/true
        mu_fine     = covguided.unitaryEspritSparse(Yb_fine, SoICols_fine, d, M);
        mu_wideFine = covguided.unitaryEspritSparse(Yb_wide, SoICols_wide, d, M);
        mu_trueFine = covguided.unitaryEspritSparse(Yb_true, SoICols_true, d, M);

        % Covariances for power estimation (fine)
        RY_fine      = (Yb_fine * Yb_fine') / N;      RY_fine      = (RY_fine      + RY_fine')/2;
        RY_wideFine  = (Yb_wide * Yb_wide') / N;      RY_wideFine  = (RY_wideFine  + RY_wideFine')/2;
        RY_trueFine  = (Yb_true * Yb_true') / N;      RY_trueFine  = (RY_trueFine  + RY_trueFine')/2;

        % power estimation
        [~, ~, mu_fineSort,     R_fine]     = covguided.estimatePowersFine(RY_fine,     SoICols_fine, mu_fine,     M);
        [~, ~, mu_wideFineSort, R_wideFine] = covguided.estimatePowersFine(RY_wideFine, SoICols_wide, mu_wideFine, M);
        [~, ~, mu_trueSort,     R_trueFine] = covguided.estimatePowersFine(RY_trueFine, SoICols_true, mu_trueFine, M);
        
        % --- build local row results (scalars -> 1x6 vectors)
        e_vec = [ ...
        sqrt(mean(angle(exp(1i*(mu_true - mu_courseSort   ))).^2)), ...
        sqrt(mean(angle(exp(1i*(mu_true - mu_wideSort     ))).^2)), ...
        sqrt(mean(angle(exp(1i*(mu_true - mu_fineSort     ))).^2)), ...
        sqrt(mean(angle(exp(1i*(mu_true - mu_wideFineSort ))).^2)), ...
        sqrt(mean(angle(exp(1i*(mu_true - mu_trueSort     ))).^2))  ...
        ];

        lpa_vec = [ ...
        covguided.computeSubspaceAngleMetrics(R_signal, R_course,    d), ...
        covguided.computeSubspaceAngleMetrics(R_signal, R_wide,      d), ...
        covguided.computeSubspaceAngleMetrics(R_signal, R_fine,      d), ...
        covguided.computeSubspaceAngleMetrics(R_signal, R_wideFine,  d), ...
        covguided.computeSubspaceAngleMetrics(R_signal, R_trueFine,  d)  ...
        ];

        % --- single consistent sliced writes
        e_trials(iter,   :) = e_vec;
        lpa_trials(iter, :) = lpa_vec;

        % update progress bar
        send(q, 1);   

    end % parfor

    %%% --- Reduce after parfor (vectorized) ---
    % RMSE^2 (MSE) per method at this ASNR
    MSE(:, snrIndex) = mean(e_trials.^2, 1).';  % since e_trials is already RMS per iter

    % Store errors/LPA for ECDF & scatter
    for m = 1:numMethods
        ERR{m, snrIndex}       = e_trials(:, m);
        LPA_TRIAL{m, snrIndex} = lpa_trials(:, m);
    end

    % Failure rates vs CRB threshold
    thr = k * max(sqrtCRB(snrIndex), 1e-12);
    failRate(:, snrIndex) = mean(e_trials > thr, 1).';

end % SNR loop

%% --- Post-processing & plots (unchanged except cosmetic) ----------------
RMSE    = sqrt(MSE);
gap_lin = RMSE ./ max(sqrtCRB, 1e-12);
gap_dB  = 20*log10(max(gap_lin, 1e-12));

%% ====== Common style ======
% Original descriptive names (for titles / text if needed)
methodNames = ["Proposed Coarse","[1] Coarse","Proposed Fine","[1] Fine","Oracle"];

% More compact legend labels (better for narrow IEEE column)
legendNames = {'Prop coarse','[1] coarse','Prop fine','[1] fine','Oracle'};

% Consistent colors (match A2/S1 style)
C = lines(numMethods);

% Markers & line styles per method
Mkr = {'o','s','^','d','>'};
LSt = {'-','-.','--',':','-'};

% Line widths, marker size, legend font size
lw    = 1.0;    % line width (match lineWidth in A2/S1)
ms    = 4;      % marker size (match markerSize)
axfs  = 8;      % axes tick labels (fontSizeAxes)
legfs = 8;      % legend font size (fontSizeLegend)

%% ====== IEEE single-column figure style (match A2/S1) ======
fontName       = 'Times New Roman';
fontSizeAxes   = 8;      % tick labels
fontSizeLabel  = 9;      % axis labels
fontSizeLegend = 7;
axisLineWidth  = 0.75;

singleColWidth = 8.8;    % cm, IEEE single-column width
figHeight      = 6.5;    % cm


%% ====== Fig 1: RMSE vs ASNR (log-y) with sqrt(CRB) ======
figure('Color','w','Units','centimeters');
set(gcf,'Position',[2 2 singleColWidth figHeight]);
set(gcf,'Renderer','painters');
hold on;
for m = 1:numMethods
    plot(ASNR, RMSE(m,:), ...
        'LineWidth', lw, ...
        'Color', C(m,:), 'LineStyle', LSt{m}, ...
        'Marker', Mkr{m}, 'MarkerSize', ms, ...
        'DisplayName', methodNames(m));
end
hCRB = plot(ASNR, sqrtCRB, '--', ...
            'Color',[0.3 0.3 0.3], ...
            'LineWidth', lw, 'DisplayName','CRB');
set(gca,'YScale','log');              % <-- force log scaling
ylim([1e-3 2]);                       % <-- choose explicit limits
yticks([1e-3 1e-2 1e-1 1]);
ytickformat('%.0e');
grid on; box on;
xlabel('ASNR [dB]', 'FontName',fontName, 'FontSize',fontSizeLabel, 'Interpreter','latex');
ylabel('RMSE [rad]', 'FontName',fontName, 'FontSize',fontSizeLabel, 'Interpreter','latex');
%
lgd_RMSE = legend(legendNames, ...
             'Location','southoutside', ...
             'Orientation','horizontal', ...
             'NumColumns', 3, ...               % <<< 3 columns → 2 rows
             'FontName', fontName, ...
             'FontSize', fontSizeLegend, ...
             'Interpreter','none', ...          % safer, narrower text
             'Box','on');
% Make legend markers/lines smaller than in the plot
legLines_RMSE = findobj(lgd_RMSE, 'Type','line');
set(legLines_RMSE, 'LineWidth', 0.75, 'MarkerSize', 3);
%
prettyAxes(gca, axfs);
exportgraphics(gcf, fullfile(figDir,'fig1_rmse_vs_asnr.pdf'), 'ContentType','vector');

%% ====== Fig 2: Gap-to-CRB vs ASNR ======
figure('Color','w','Units','centimeters');
set(gcf,'Position',[2 2 singleColWidth figHeight]);
set(gcf,'Renderer','painters');
hold on;
for m = 1:numMethods
    plot(ASNR, gap_dB(m,:), 'LineWidth', lw, ...
        'Color', C(m,:), 'LineStyle', LSt{m}, 'Marker', Mkr{m}, 'MarkerSize', ms, ...
        'DisplayName', methodNames(m));
end
%
% --- draw unlabeled reference line ---
yline(0, '--', ...
      'Color',[0.4 0.4 0.4], ...
      'LineWidth', 1.0, ...
      'HandleVisibility','off');

% --- add manual label near left side (inside axes) ---
ax = gca;
xl = xlim(ax);

xText = xl(1) + 0.03*(xl(2) - xl(1));   % 3% in from the left
yText = 0;                              % same height as the line


text(xText, yText, 'sqrt(CRB)', ...
     'Parent', ax, ...
     'HorizontalAlignment','left', ...
     'VerticalAlignment','bottom', ...
     'FontName', fontName, ...
     'FontSize', fontSizeLabel, ...
     'Interpreter','none');   % or 'tex'

%
grid on; box on;
xlabel('ASNR [dB]', 'FontName',fontName, 'FontSize',fontSizeLabel, 'Interpreter','latex');
ylabel('Gap to CRB [dB]', 'FontName',fontName, 'FontSize',fontSizeLabel, 'Interpreter','latex');
lgd_GapCRB = legend(legendNames, ...
             'Location','southoutside', ...
             'Orientation','horizontal', ...
             'NumColumns', 3, ...               % <<< 3 columns → 2 rows
             'FontName', fontName, ...
             'FontSize', fontSizeLegend, ...
             'Interpreter','none', ...          % safer, narrower text
             'Box','on');
% Make legend markers/lines smaller than in the plot
legLines_GapCRB = findobj(lgd_GapCRB, 'Type','line');
set(legLines_GapCRB, 'LineWidth', 0.75, 'MarkerSize', 3);

prettyAxes(gca, axfs);
exportgraphics(gcf, fullfile(figDir,'fig2_gap_to_crb.pdf'), 'ContentType','vector');

%% ====== Fig 3: Failure-Rate vs ASNR with Wilson 95% CI (linear scale, % units) ======
figure('Color','w','Units','centimeters');
set(gcf,'Position',[2 2 singleColWidth figHeight]);
set(gcf,'Renderer','painters');
hold on;
%
n = totalIter;                           % trials per ASNR
scaleVis = 1;                            % =1 exact; >1 widens bands for visibility

lineHandles = gobjects(1,numMethods);
for m = 1:numMethods
    p  = failRate(m,:);                       % fraction in [0,1]
    [lo, hi] = wilsonCI(p, n, 0.95);          % CI in [0,1]

    % --- optional visual widening (purely cosmetic, doesn't touch mean curve)
    if scaleVis ~= 1
      lo = max(0, p - scaleVis*(p - lo));
      hi = min(1, p + scaleVis*(hi - p));
    end

    x = ASNR;

    % --- shaded band
    hBand = patch('XData',[x fliplr(x)], ...
                  'YData',100*[lo fliplr(hi)], ...
                  'FaceColor', C(m,:), 'FaceAlpha', 0.18, ...
                  'EdgeColor','none', 'HandleVisibility','off');

    % --- dashed CI bounds (make them clearly visible)
    hLo = plot(x, 100*lo, '--', 'Color', C(m,:)*0.55, 'LineWidth', 1.2, 'HandleVisibility','off');
    hHi = plot(x, 100*hi, '--', 'Color', C(m,:)*0.55, 'LineWidth', 1.2, 'HandleVisibility','off');

    % --- mean curve on top
    lineHandles(m) = plot(x, 100*p, ...
        'LineWidth', lw, 'Color', C(m,:), ...
        'LineStyle', LSt{m}, 'Marker', Mkr{m}, 'MarkerSize', ms, ...
        'DisplayName', methodNames(m));

    % ensure mean line is above band/edges
    uistack(lineHandles(m), 'top');
end

% --- axes formatting (log in percent)
ylim([0 100]);                              % 0% to 100%
yticks(0:20:100);
xlabel('ASNR [dB]', 'FontName',fontName, 'FontSize',fontSizeLabel, 'Interpreter','latex');
ylabel('Failure rate [\%]', 'FontName',fontName, 'FontSize',fontSizeLabel, 'Interpreter','latex');
grid on; box on;
lgd_Fail = legend(lineHandles, legendNames, ...
             'Location','southoutside', ...
             'Orientation','horizontal', ...
             'NumColumns', 3, ...
             'FontName', fontName, ...
             'FontSize', fontSizeLegend, ...
             'Interpreter','none', ...
             'Box','on');
% Make legend markers/lines smaller than in the plot
legLines_Fail = findobj(lgd_Fail, 'Type','line');
set(legLines_Fail, 'LineWidth', 0.75, 'MarkerSize', 3);

prettyAxes(gca, axfs);

% optional: vector export
exportgraphics(gcf, fullfile(figDir,'fig3_failrate.pdf'), 'ContentType','vector');

%% ====== Fig 4: ECDFs at 0 dB and 15 dB (Proposed Fine, [5] Fine, Oracle) ======
lowIdx = find(ASNR == -5, 1);
midIdx = find(ASNR == 0, 1);
midHighIdx = find(ASNR == 5, 1);
hiIdx  = find(ASNR == 15, 1);
snrIdxs = [lowIdx, midIdx, midHighIdx, hiIdx];
labelsECDF = {'Proposed Fine','[1] Fine','Oracle'};
methECDF   = [3, 4, 5];  % indices in your ERR cell

%
figure('Color','w','Units','centimeters');
set(gcf,'Position',[2 2 singleColWidth figHeight]);
set(gcf,'Renderer','painters');

tl = tiledlayout(2,2, 'Padding','compact', 'TileSpacing','compact');

%
hECDF = gobjects(1, numel(methECDF));   % handles for shared legend

for t = 1:length(snrIdxs)
    s = snrIdxs(t);
    ax = nexttile;                      % creates the next subplot
    hold(ax, 'on');

    for j = 1:numel(methECDF)
        m = methECDF(j);
        [f, x] = ecdf(ERR{m, s});

        if t == 1
            % Store handles from the FIRST subplot for the global legend
            hECDF(j) = plot(ax, x, f, 'LineWidth', lw, ...
                'Color', C(m,:), 'LineStyle', LSt{m}, ...
                'DisplayName', labelsECDF{j});
        else
            % Same curves in other subplots but hidden from legend
            plot(ax, x, f, 'LineWidth', lw, ...
                'Color', C(m,:), 'LineStyle', LSt{m}, ...
                'HandleVisibility', 'off');
        end
    end

    grid(ax, 'on'); box(ax, 'on');
    if (t==3)
        xlabel(ax, 'Per-trial RMSE [rad]', ...
           'FontName',fontName, 'FontSize',fontSizeLabel, 'Interpreter','latex');
    end
    if (t==1) || (t==3) 
        ylabel(ax, 'ECDF', ...
           'FontName',fontName, 'FontSize',fontSizeLabel, 'Interpreter','latex');
    end
    title(ax, sprintf('ASNR = %.0f dB', ASNR(s)), ...
          'FontName',fontName, 'FontSize',fontSizeLabel, 'Interpreter','latex');
    if t==3
        xlim(ax, [-0.01 0.2]);
    end
    prettyAxes(ax, axfs);
end
% legend
lgd_ECDF = legend(hECDF, labelsECDF, ...
             'Orientation','horizontal', ...
             'NumColumns', numel(labelsECDF), ...
             'FontName', fontName, ...
             'FontSize', fontSizeLegend, ...
             'Interpreter','latex', ...
             'Box','on');

% Attach the legend as its own tile below the 2×2 grid
lgd_ECDF.Layout.Tile = 'south';

% export graphics
exportgraphics(gcf, fullfile(figDir,'fig4_ecdf_0dB_15dB.pdf'), 'ContentType','vector');

%% ====== Fig 5: LPA–Error Scatter at 15 dB (Proposed Fine vs [5] Fine) ======
snrScatter = hiIdx;   % 15 dB
figure('Color','w','Units','centimeters');
set(gcf,'Position',[2 2 singleColWidth figHeight]);
set(gcf,'Renderer','painters');
hold on;
mList = [3, 4];  % Proposed Fine, [5] Fine
for ii = 1:numel(mList)
    m = mList(ii);
    e  = ERR{m, snrScatter};       % per-trial RMS error [rad]
    a  = LPA_TRIAL{m, snrScatter}; % per-trial LPA [deg]
    scatter(a, e, 12, 'MarkerFaceColor', C(m,:), 'MarkerEdgeColor', C(m,:), ...
        'MarkerFaceAlpha', 0.15, 'MarkerEdgeAlpha', 0.15, ...
        'DisplayName', methodNames(m));
end
xlabel('Largest Principal Angle [deg]', 'FontName',fontName, 'FontSize', fontSizeLabel, 'Interpreter','latex');
ylabel('Per-trial RMSE [rad]', 'FontName',fontName, 'FontSize', fontSizeLabel, 'Interpreter','latex');
legend('Location','northeast', ...
       'FontName',fontName, ...
       'FontSize',fontSizeLegend, ...
       'Interpreter','latex', ...
       'Box','on');
grid on; box on;
prettyAxes(gca, axfs);
exportgraphics(gcf, fullfile(figDir,'fig5_scatter_15dB.pdf'), 'ContentType','vector');

%% ====== Local helpers (keep at end of script) ======
function [lo, hi] = wilsonCI(p, n, conf)
% p: row vector of proportions (0..1), n: scalar trials, conf: e.g., 0.95
% returns lo, hi (row vectors) for Wilson score interval
if nargin < 3, conf = 0.95; end
z = -sqrt(2)*erfcinv(conf);  % inverse of N(0,1) CDF for two-sided conf
p = p(:).';                  % force row
den = 1 + (z^2)/n;
center = (p + (z^2)/(2*n)) ./ den;
half = (z ./ den) .* sqrt( (p.*(1-p)/n) + (z^2)/(4*n^2) );
lo = max(0, center - half);
hi = min(1, center + half);
end

function prettyAxes(ax, fs)
% Consistent axis cosmetics (IEEE single-column style)
set(ax, 'LineWidth', 0.75, ...              % axisLineWidth
        'FontSize', fs, ...                 % use axfs = 8
        'FontName', 'Times New Roman', ...
        'TickDir','out', 'Box','on', ...
        'MinorGridLineStyle','-', ...
        'GridAlpha', 0.2, ...
        'TickLabelInterpreter','latex');
ax.XGrid = 'on';
ax.YGrid = 'on';
end


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
%% ---------- Aggregates ----------
% method_names = ["Proposed Coarse","[5] Coarse","Proposed Fine","[5] Fine","Oracle"];
% [MM, SS] = size(MSE);
% gap_dB = 20*log10(max(MSE.^0.5 ./ max(sqrtCRB,1e-12), 1e-12)); % optional
% 
% Agg = table;
% for m = 1:MM
%     Agg = [Agg; table( ...
%         ASNR(:), repmat(m,SS,1), repmat(method_names(m),SS,1), ...
%         sqrt(MSE(m,:)).', failRate(m,:).', sqrtCRB(:), gap_dB(m,:).', ...
%         'VariableNames', {'ASNR','method','method_name','RMSE','FailRate','sqrtCRB','gap_dB'})]; %#ok<AGROW>
% end
% writetable(Agg, fullfile(csvDir,'aggregates.csv'));
%% ---------- Per-trial subset for ECDF/scatter ----------
% wantASNR = [-5, 0, 5, 15];           % adjust if you like
% Trials = table;
% for a = wantASNR
%     s = find(ASNR == a, 1);
%     if isempty(s), warning('ASNR %g dB not found, skipping', a); continue; end
%     for m = 1:MM
%         e = ERR{m,s};           % vector totalIter x 1 (per-trial RMS error)
%         l = LPA_TRIAL{m,s};     % vector totalIter x 1 (LPA in degrees)
%         Trials = [Trials; table( ...
%             repmat(a, numel(e),1), repmat(m, numel(e),1), (1:numel(e)).', ...
%             e(:), l(:), repmat(method_names(m), numel(e),1), ...
%             'VariableNames', {'ASNR','method','trial','err_rms','lpa_deg','method_name'})]; %#ok<AGROW>
%     end
% end
% writetable(Trials, fullfile(csvDir,'trials.csv'));