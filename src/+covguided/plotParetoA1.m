function plotParetoA1(med_t_total_ms, rmse_rad, varargin)
%PLOT_PARETO_A1  Pareto RMSE vs Runtime (A1) with per-SNR hulls + optional dashed guides + side table.
% Styled to match A2 figures (Times New Roman, IEEE single-column, lines(3) palette).
%
% Usage:
%   plot_pareto_A1(med_t_total_ms, rmse_rad, ...
%       'snr', [3 6], ...
%       'condLabels', {'Cov 12->6','Cov 12->12','Sect 12->6','Sect 12->12'}, ...
%       'breakdown', bk, ...
%       'connect_all_per_snr', true, ...
%       'print_hull', true, ...
%       'outpdf', fullfile('figures_pdf','A1_pareto.pdf'), ...
%       'tablecsv', fullfile('figures_pdf','A1_timing_table.csv'));

% ---------- parse inputs ----------
parse = inputParser;
parse.addParameter('snr',[3 6], @(x)isnumeric(x)&&isvector(x));
parse.addParameter('condLabels',{}, @(x)iscellstr(x) || isstring(x));
parse.addParameter('style', struct(), @(x)isstruct(x));
parse.addParameter('breakdown', struct(), @(x)isstruct(x));
parse.addParameter('outpdf', 'A1_pareto.pdf', @(x)ischar(x) || isstring(x));
parse.addParameter('tablecsv', 'A1_timing_table.csv', @(x)ischar(x) || isstring(x));
parse.addParameter('connect_all_per_snr', false, @(b)islogical(b)&&isscalar(b));
parse.addParameter('print_hull', true, @(b)islogical(b)&&isscalar(b));
parse.parse(varargin{:});

snrList   = parse.Results.snr;
condLbl   = parse.Results.condLabels;
sty       = parse.Results.style;
bk        = parse.Results.breakdown;
outpdf    = char(parse.Results.outpdf);
tablecsv  = char(parse.Results.tablecsv);
doAll     = parse.Results.connect_all_per_snr;
doPrint   = parse.Results.print_hull;

[C,S] = size(med_t_total_ms);
assert(isequal(size(rmse_rad), [C,S]), 'rmse_rad must be CxS matching med_t_total_ms.');
if isempty(condLbl), condLbl = compose("Cond %d", 1:C); end

% ---------- visual style (match A2) ----------
fontName       = 'Times New Roman';
fontSizeAxes   = 8;      % tick labels
fontSizeLabel  = 9;      % axis labels
fontSizeLegend = 8;
lineWidth      = 1.0;
axisLineWidth  = 0.75;
markerSize     = 4;

singleColWidth = 8.8;    % cm, IEEE single column
figHeight      = 6.5;    % cm

colors = lines(3);       % 1: cov-guided, 2: sector, 3: spare/CRB

% Allow optional override via sty, but default to lines(3)
if ~isfield(sty,'covColor'),    sty.covColor    = colors(1,:); end
if ~isfield(sty,'sectColor'),   sty.sectColor   = colors(2,:); end
if ~isfield(sty,'shape12to6'),  sty.shape12to6  = 'o'; end
if ~isfield(sty,'shape12to12'), sty.shape12to12 = 's'; end

% ---------- map conditions to color / marker ----------
isCov  = contains(string(condLbl),'Cov','IgnoreCase',true) | contains(string(condLbl),'Guided','IgnoreCase',true);
is126  = contains(string(condLbl),'12->6')  | contains(string(condLbl),'12\u21926');
is1212 = contains(string(condLbl),'12->12') | contains(string(condLbl),'12\u219212');

col = zeros(C,3);
mk  = repmat('o',C,1);
for c=1:C
    if isCov(c), col(c,:) = sty.covColor; else, col(c,:) = sty.sectColor; end
    if is126(c)
        mk(c) = sty.shape12to6;
    elseif is1212(c)
        mk(c) = sty.shape12to12;
    else
        mk(c) = 'o';
    end
end

% ---------- figure / axes setup ----------
fig = figure('Color','w', ...
             'Units','centimeters', ...
             'Visible','off');
fig.Position(3:4) = [singleColWidth, figHeight];
set(fig, 'Renderer','painters');

ax = axes(fig);
hold(ax, 'on'); grid(ax, 'on'); box(ax, 'on');

% ---------- plot points (markers only) ----------
hPts = gobjects(C,1);
for c=1:C
    x = med_t_total_ms(c,:);
    y = rmse_rad(c,:);
    hPts(c) = plot(ax, x, y, ...
        'LineStyle','none', ...
        'Marker', mk(c), ...
        'MarkerSize', markerSize, ...
        'LineWidth', lineWidth, ...
        'MarkerFaceColor', col(c,:), ...
        'MarkerEdgeColor', [0 0 0], ...
        'Color', col(c,:));  % Color used if we ever add per-condition lines
end

% ---------- optional dashed connectors per SNR ----------
if doAll
    for s=1:S
        [xs,ord] = sort(med_t_total_ms(:,s), 'ascend');
        ys = rmse_rad(ord,s);
        plot(ax, xs, ys, '--', ...
            'LineWidth', lineWidth*0.8, ...
            'Color', [0.7 0.7 0.7], ...
            'HandleVisibility','off');
    end
end

%% ---------- Pareto hull per SNR (strict) ----------
hHull = gobjects(S,1);  % one handle per SNR hull

for s=1:S
    x = med_t_total_ms(:,s);
    y = rmse_rad(:,s);
    idx = pareto_front(x,y);  % local helper below

    if doPrint
        fprintf('Pareto @ %gdB -> rows: %s\n', snrList(min(s,end)), mat2str(find(idx).'));
    end

    [xs,ord] = sort(x(idx)); ys = y(idx); ys = ys(ord);

    % Style per SNR: inner = dark gray (3 dB), outer = black (6 dB)
    if s == 1
        hullColor = [0.3 0.3 0.3];  % 3 dB
    else
        hullColor = [0 0 0];        % 6 dB
    end

    hHull(s) = plot(ax, xs, ys, '-', ...
        'LineWidth', lineWidth*1.5, ...
        'Color', hullColor);
end

%% ---------- axes labels / style ----------
xlabel(ax, 'Runtime per trial [ms]', ...
       'FontName', fontName, ...
       'FontSize', fontSizeLabel, ...
       'Interpreter','latex');
ylabel(ax, 'RMSE [rad]', ...
       'FontName', fontName, ...
       'FontSize', fontSizeLabel, ...
       'Interpreter','latex');

set(ax, 'FontName', fontName, ...
        'FontSize', fontSizeAxes, ...
        'LineWidth', axisLineWidth, ...
        'TickLabelInterpreter','latex');

xlim(ax, [0, max(med_t_total_ms(:))*1.1]);
ylim(ax, [0, max(rmse_rad(:))*1.1]);

%% ---------- legend (single, horizontal, below) ----------
% Build hull labels explicitly, e.g. "3 dB hull", "6 dB hull"
% hullLabels = arrayfun(@(x) sprintf('%g dB hull', x), snrList(:), 'UniformOutput', false);
% 
% % Stack condition markers first, then hull lines
% legHandles = [hPts(:); hHull(:)];
% legLabels  = [condLbl(:); hullLabels(:)];
% 
% lgd = legend(ax, legHandles, legLabels, ...
%     'Orientation','horizontal', ...
%     'NumColumns', 2, ...              % 2 columns: 3 rows x 2 columns (6 entries)
%     'Box','on', ...
%     'FontName', fontName, ...
%     'FontSize', fontSizeLegend, ...
%     'Interpreter','none', ...         % show "->" literally
%     'Location','southoutside');
% 
% lgd.ItemTokenSize = [10 6];

% ---------- legend (inside, top-left, compact) ----------
hullLabels = arrayfun(@(x) sprintf('%g dB hull', x), snrList(:), 'UniformOutput', false);

legHandles = [hPts(:); hHull(:)];
legLabels  = [condLbl(:); hullLabels(:)];

lgd = legend(ax, legHandles, legLabels, ...
    'Orientation','vertical', ...
    'NumColumns', 1, ...                % single column, narrow
    'Box','on', ...
    'FontName', fontName, ...
    'FontSize', fontSizeLegend-1, ...   % slightly smaller
    'Interpreter','none', ...
    'Location','northwest');            % inside, top-left

lgd.ItemTokenSize = [8 6];

%% ---------- save vector PDF ----------
if ~isempty(outpdf)
    try
        if exist('covguided.saveVectorPdf','file')==2 || exist('saveVectorPdf','file')==2
            covguided.saveVectorPdf(fig, outpdf);
        else
            exportgraphics(fig, outpdf, 'ContentType','vector','BackgroundColor','none');
        end
        fprintf('Saved vector PDF to %s\n', outpdf);
    catch ME
        warning('plot_pareto_A1:SaveFailed','[%s] %s', ME.identifier, ME.message);
    end
end

%% ---------- side timing table ----------
Cond     = strings(0,1);
SNR      = zeros(0,1);
RMSE_col = [];        % NEW: RMSE per condition/SNR
t_cov    = [];
t_sel    = [];
t_es     = [];
t_tot    = [];
iters    = [];
rfct     = [];

hasBk = isfield(bk,'cov_ms') && isfield(bk,'sel_ms') && ...
        isfield(bk,'es_ms')  && isfield(bk,'total_ms');

if hasBk
    for c = 1:C
        for s = 1:S
            Cond(end+1,1)     = string(condLbl{c});
            SNR(end+1,1)      = snrList(min(s, numel(snrList)));
            RMSE_col(end+1,1) = rmse_rad(c,s);       % <-- NEW

            t_cov(end+1,1)    = bk.cov_ms(c,s);
            t_sel(end+1,1)    = bk.sel_ms(c,s);
            t_es(end+1,1)     = bk.es_ms(c,s);
            t_tot(end+1,1)    = bk.total_ms(c,s);

            if isfield(bk,'iters')
                iters(end+1,1) = bk.iters(c,s);
            else
                iters(end+1,1) = NaN;
            end
            if isfield(bk,'rf_chain_time')
                rfct(end+1,1) = bk.rf_chain_time(c,s);
            else
                rfct(end+1,1) = NaN;
            end
        end
    end

    T = table(Cond, SNR, RMSE_col, t_cov, t_sel, t_es, t_tot, iters, rfct, ...
        'VariableNames', {'Condition','SNR_dB','RMSE_rad', ...
                          't_cov_ms','t_sel_ms','t_es_ms','t_total_ms', ...
                          'iters','rf_chain_time'});

    disp('--- Median timing/accuracy table (per condition x SNR) ---');
    disp(T);
    if ~isempty(tablecsv)
        try
            writetable(T, tablecsv);
            fprintf('Saved timing table CSV to %s\n', tablecsv);
        catch ME
            warning('plot_pareto_A1:TableWriteFailed','[%s] %s', ...
                    ME.identifier, ME.message);
        end
    end
else
    disp('No breakdown struct provided; skipping timing table CSV.');
end


close(fig);  % mirror A2 "Visible','off' + export + close"

end % main function


%% ---------- helpers ----------
function idx = pareto_front(x, y)
%PARETO_FRONT Return logical index of points on the non-dominated front.
% Minimize both objectives x and y.
n = numel(x);
idx = false(n,1);
[~, ord] = sort(x, 'ascend');
besty = inf;
for k = 1:n
    i = ord(k);
    if y(i) < besty - eps
        idx(i) = true;
        besty = y(i);
    end
end
end