function runPhaseQuantizationStudy(varargin)
%RUNPHASEQUANTIZATIONSTUDY  R1.1 Monte Carlo: graceful degradation under
%   B-bit phase-shifter quantization.
%
%   Mirrors runMainExperiment.m exactly, with one physical modification:
%   the fine-stage analog combiner is a B-bit quantized DFT dictionary
%   Bq in place of the ideal DFT. The coarse stage uses element-space
%   subsampling (no phase shifters) and is invariant to B, so all
%   observed RMSE changes are attributable to fine-stage quantization.
%
%   The ideal-DFT fast path `phase_k .* fft(Y,M,1)` is replaced by the
%   explicit product `Bq' * Y_fine` for each B. Honesty rule: every B
%   in BitsList is written to CSV, no cherry-picking.
%
%   Usage:
%     runPhaseQuantizationStudy();            % defaults
%     runPhaseQuantizationStudy('Trials',20,'ASNRdB',[0 15],'BitsList',[2 4 Inf]);

addpath(genpath(fullfile(fileparts(mfilename('fullpath')), '..', 'src')));
addpath(fullfile(fileparts(mfilename('fullpath')), '..', 'config'));

p = inputParser;
addParameter(p,'Trials',1e4,@(x)isscalar(x)&&x>0);
addParameter(p,'ASNRdB',[0 5 10 15]);
addParameter(p,'BitsList',[1 2 3 4 5 6 Inf]);
parse(p,varargin{:});
opt = p.Results;

% ---- invariants from runMainExperiment.m --------------------------------
sim = getSimParams();  %#ok<NASGU>   % ensures config path resolves
M   = 32;  N = 100;  NRF = 12;  d = 3;
pows     = [0.95, 0.5, 0.1];
R_source = diag(pows);
Psqrt    = diag(sqrt(pows/2));
mu_true  = [-2.1; 0.5; 2.5];
Atrue    = exp(1i * (0:M-1).' * mu_true.');
R_signal = Atrue * R_source * Atrue';

maskIdx  = covguided.centeredContiguousMask(M, NRF);
mIdxFlip = NRF:-1:1;

Bideal   = covguided.dftDictionary(M);

repoRoot = fullfile(fileparts(mfilename('fullpath')), '..');
figDir = fullfile(repoRoot,'results','figures');
csvDir = fullfile(repoRoot,'results','csv');
if ~exist(figDir,'dir'), mkdir(figDir); end
if ~exist(csvDir,'dir'), mkdir(csvDir); end

asnrList = opt.ASNRdB;
bitsList = opt.BitsList(:).';
noiseVar = real(trace(R_signal))./(M*10.^(asnrList/10));
sqrtCRB  = localSqrtCRB(Atrue, R_source, noiseVar, N, d);

Rt = opt.Trials;  nB = numel(bitsList);  nS = numel(asnrList);
kThr = 3;

% Precompute quantized dictionaries
BqCell = cell(1,nB);
for iB = 1:nB, BqCell{iB} = quantizeDictionary(Bideal, bitsList(iB)); end

MSE      = zeros(nB,nS);
FailRate = zeros(nB,nS);

fprintf('[R1.1] Phase-quant study: Trials=%d, B-values=%d, ASNRs=%d\n', ...
        Rt, nB, nS);

for iS = 1:nS
    nvStd = sqrt(noiseVar(iS)/2);
    fThr  = kThr * sqrtCRB(iS);

    e_trials    = zeros(Rt,nB);
    fail_trials = false(Rt,nB);

    parfor iter = 1:Rt
        s = RandStream('Threefry','Seed',5489);
        s.Substream = iter;
        RandStream.setGlobalStream(s);

        % --- Coarse stage (B-invariant, element-space) ---------------
        S  = Psqrt * (randn(d,N) + 1i*randn(d,N));
        Z  = nvStd * (randn(M,N) + 1i*randn(M,N));
        Y_coarse  = Atrue*S + Z;

        Yb_coarse = Y_coarse(maskIdx, :);
        RY_coarse = (Yb_coarse*Yb_coarse')/N;
        RY_coarse = (RY_coarse+RY_coarse')/2;
        RY_fba    = (RY_coarse + conj(RY_coarse(mIdxFlip, mIdxFlip))) * 0.5;

        mu_coarse = covguided.tlsEsprit(RY_fba, d);
        [~,~,mu_courseSort, R_course] = covguided.estimatePowersCoarse( ...
                                          RY_fba, maskIdx, mu_coarse, M);

        beamGroups_fine = covguided.sectorizeDftBeams(mu_courseSort, M, 4, 1);
        selBeams_fine   = covguided.selectAdjacentPairsFromSectorizedCov( ...
                                          R_course, beamGroups_fine);
        SoICols_fine    = unique(cell2mat(selBeams_fine)).';

        % --- Fine-stage data (B-invariant draws) ---------------------
        S2 = Psqrt * (randn(d,N) + 1i*randn(d,N));
        Z2 = nvStd * (randn(M,N) + 1i*randn(M,N));
        Y_fine = Atrue*S2 + Z2;

        % --- Inner B-loop: quantization enters ONLY here -------------
        e_row = zeros(1,nB);  f_row = false(1,nB);
        for iB = 1:nB
            Bq = BqCell{iB};
            Yb_full = Bq' * Y_fine;
            Yb_fine = Yb_full(SoICols_fine, :);

            mu_fine = covguided.unitaryEspritSparse( ...
                          Yb_fine, SoICols_fine, d, M);

            RY_fine = (Yb_fine*Yb_fine')/N;
            RY_fine = (RY_fine+RY_fine')/2;
            [~,~,mu_fineSort,~] = covguided.estimatePowersFine( ...
                                      RY_fine, SoICols_fine, mu_fine, M);

            err = angle(exp(1i*(mu_true - mu_fineSort)));
            e_row(iB) = sqrt(mean(err.^2));
            f_row(iB) = e_row(iB) > fThr;
        end
        e_trials(iter,:)    = e_row;
        fail_trials(iter,:) = f_row;
    end

    MSE(:,iS)      = mean(e_trials.^2, 1).';
    FailRate(:,iS) = mean(fail_trials, 1).';
    fprintf('  ASNR=%2d dB: RMSE(B=Inf)=%.3e, RMSE(B=%g)=%.3e\n', ...
            asnrList(iS), sqrt(MSE(end,iS)), bitsList(1), sqrt(MSE(1,iS)));
end

RMSE  = sqrt(MSE);
bound = arrayfun(@(B) sin(pi/2^max(B,1)), bitsList);
bound(isinf(bitsList)) = 0;

% ---- CSV ---------------------------------------------------------------
rows = {};
for iB = 1:nB
    Btag = bitsList(iB);
    Bstr = sprintf('Cov_B%s', ternary(isinf(Btag),'Inf',num2str(Btag)));
    Benc = Btag;  if isinf(Benc), Benc = -1; end
    for iS = 1:nS
        rows(end+1,:) = {2, asnrList(iS), Bstr, RMSE(iB,iS), ...
                         FailRate(iB,iS), sqrtCRB(iS), kThr, ...
                         Rt, Benc, bound(iB)}; %#ok<AGROW>
    end
end
T = cell2table(rows, 'VariableNames', {'Kf','ASNR_dB','method', ...
    'RMSE_rad','FailRate','sqrtCRB_rd','kThr','totalIter', ...
    'Bbits','AnalyticalBound_l2'});

stamp   = datestr(now,'yyyymmdd_HHMMSS');
csvPath = fullfile(csvDir, sprintf('phase_quant_study_%s.csv', stamp));
writetable(T, csvPath);
fprintf('[R1.1] wrote %s (%d rows)\n', csvPath, height(T));

% ---- Figure ------------------------------------------------------------
fig = figure('Visible','off','Units','centimeters','Position',[0 0 9 9]);  % taller
tl  = tiledlayout(fig, 1, 1, 'Padding','compact', 'TileSpacing','compact');
ax  = nexttile(tl); hold(ax,'on'); set(ax,'YScale','log'); grid(ax,'on');
cmap = lines(nB);
for iB = 1:nB
    lbl = sprintf('B=%s', ternary(isinf(bitsList(iB)),'\infty', ...
                                  num2str(bitsList(iB))));
    plot(ax, asnrList, RMSE(iB,:), '-o', 'Color',cmap(iB,:), ...
         'LineWidth',1.1,'DisplayName',lbl);
end
plot(ax, asnrList, sqrtCRB, 'k--','LineWidth',1.0,'DisplayName','\surdCRB');
xlabel(ax,'ASNR [dB]'); ylabel(ax,'RMSE [rad]');
title(ax,'Phase-quantization degradation (R1.1)');
lgd = legend(ax, 'Orientation','horizontal','NumColumns',4);
lgd.Layout.Tile = 'south';            % pushes legend below axes
hold(ax,'off');
figPath = fullfile(figDir, sprintf('phase_quant_rmse_%s.pdf', stamp));
covguided.saveVectorPdf(fig, figPath);  close(fig);
fprintf('[R1.1] wrote %s\n', figPath);

end

% =======================================================================
function Bq = quantizeDictionary(B, bits)
    if isinf(bits), Bq = B; return; end
    L    = 2^bits;
    phiQ = (2*pi/L) * round(L*angle(B)/(2*pi));
    Bq   = abs(B) .* exp(1i*phiQ);
end

function c = localSqrtCRB(Atrue, R_source, noiseVar, N, d)
    M = size(Atrue,1);
    PorthA = eye(M) - Atrue*pinv(Atrue);
    derA   = (1i*transpose(0:M-1)) .* Atrue;
    H      = transpose(derA' * PorthA * derA);
    CRB    = zeros(1,numel(noiseVar));
    for k = 1:numel(noiseVar)
        secondTerm = ((Atrue'*Atrue) * R_source) / noiseVar(k);
        firstTerm  = eye(d) + secondTerm;
        CRBinv     = diag(inv(real((R_source * (firstTerm\secondTerm)) .* H)));
        CRB(k)     = sum((noiseVar(k)/(2*N)) * CRBinv) / d;
    end
    c = sqrt(CRB);
end

function y = ternary(c,a,b)
    if c, y = a; else, y = b; end
end
