function [selected_pairs, meta] = select_adjacent_pairs_from_sectorized_cov(R_model, beam_groups, opts)
% Select best adjacent DFT beam pair per sector from a reconstructed covariance.
% Simplifications:
%   - L is fixed to 2 (adjacent beams only).
%   - No cross-sector duplicate/overlap restrictions.
%
% Inputs
%   R_model     : [M x M] Hermitian PSD covariance (e.g., A*Rs*A')
%   beam_groups : 1xG cell, each a vector of beam indices (sectorize_dft_beams output)
%   opts        : (optional) struct with fields:
%       .eig_tol_scale   (1e-6)   threshold for signal-rank detection
%       .alpha           (0.5)    base weight for cond^2 penalty
%       .q               (0.90)   quantile for adaptive penalty
%       .lambda_bounds   [1e-4,10]
%       .tikhonov_gamma  (0)      projector reg: (B'B + gamma I)^(-1)
%       .normalize_score (true)   divide preservation by d
%
% Outputs
%   selected_pairs : 1xG cell, each a 1x2 vector [b1 b2] (global beam indices)
%   meta           : struct with diagnostics per sector

    if nargin < 3 || isempty(opts), opts = struct(); end
    if ~isfield(opts,'eig_tol_scale')   || isempty(opts.eig_tol_scale),   opts.eig_tol_scale   = 1e-6;  end
    if ~isfield(opts,'alpha')           || isempty(opts.alpha),           opts.alpha           = 0.5;   end
    if ~isfield(opts,'q')               || isempty(opts.q),               opts.q               = 0.90;  end
    if ~isfield(opts,'lambda_bounds')   || isempty(opts.lambda_bounds),   opts.lambda_bounds   = [1e-4, 10]; end
    if ~isfield(opts,'tikhonov_gamma')  || isempty(opts.tikhonov_gamma),  opts.tikhonov_gamma  = 0;     end
    if ~isfield(opts,'normalize_score') || isempty(opts.normalize_score), opts.normalize_score = true;  end

    G = numel(beam_groups);
    M = size(R_model,1);
    selected_pairs = cell(G,1);
    meta = struct('lambda', nan(G,1), 'best_score', nan(G,1), 'best_kappa', nan(G,1));

    % Phase-shifted DFT columns (1-based bins -> 0-based phase)
    Bcols = @(bins0) exp(1i * 2*pi * (0:M-1).' .* (bins0 / M)) ...
                   .* exp(-1i * (M-1) * pi * (bins0 / M));

    % --- Denoised signal subspace (soft-threshold proxy)
    R = (R_model + R_model')/2;
    [U_all, S_vals] = eig(R, 'vector');
    [S_sorted, idx] = sort(S_vals, 'descend'); U_all = U_all(:,idx);

    tol = opts.eig_tol_scale * norm(R, 'fro');
    d = sum(S_sorted > max(tol,0));
    if d < 1, error('Detected signal rank is zero.'); end

    tau = 0; if d < numel(S_sorted), tau = median(S_sorted(d+1:end)); end
    lam = S_sorted(1:d);
    shrink = max(lam - tau, 0);
    sc = ones(d,1);
    nz = (lam > 0) & (shrink > 0);
    sc(nz) = sqrt(shrink(nz)./lam(nz));
    U_s = U_all(:,1:d) * diag(sc);
    norm_factor = d; if ~opts.normalize_score, norm_factor = 1; end

    % --- Per-sector: score all adjacent pairs; pick the best (no dedup across sectors)
    for g = 1:G
        idxs = sort(beam_groups{g}(:));
        nb = numel(idxs);

        if nb < 2
            warning('Sector %d has <2 beams; returning what is available.', g);
            selected_pairs{g} = idxs(:).';
            continue;
        end

        % Circular adjacent pairs within the sector
        wrap = [idxs; idxs(1)];     % pairs: (b1,b2),...,(b_nb,b1)
        num_pairs = nb;

        cond_vals = nan(num_pairs,1);
        pres_vals = nan(num_pairs,1);

        % First pass: kappa and preservation for each adjacent pair
        for w = 1:num_pairs
            pair = sort([wrap(w), wrap(w+1)]);   % 1x2 global indices
            Bs   = Bcols(pair-1);
            Gram = real(Bs' * Bs);

            kappaG = cond(Gram);
            if ~isfinite(kappaG) || kappaG > 1e12
                cond_vals(w) = Inf; pres_vals(w) = -Inf; continue;
            end
            cond_vals(w) = sqrt(max(kappaG,1)); % approx cond(B) via sqrt(cond(Gram))

            if opts.tikhonov_gamma > 0
                P = Bs / (Gram + opts.tikhonov_gamma*eye(2)) * Bs';
            else
                P = Bs / Gram * Bs';
            end
            pres_vals(w) = real(trace(U_s' * P * U_s)) / norm_factor; % ~[0,1] if normalized
        end

        finite_mask = isfinite(cond_vals) & isfinite(pres_vals);
        if ~any(finite_mask)
            % Fallback: per-beam energy, pick two strongest
            e = zeros(nb,1);
            for i=1:nb
                v = Bcols(idxs(i)-1);
                e(i) = real(v' * R * v);
            end
            [~,o] = sort(e,'descend');
            selected_pairs{g} = sort(idxs(o(1:min(2,nb)))).';
            meta.best_score(g) = -inf; meta.best_kappa(g) = inf; meta.lambda(g) = NaN;
            continue;
        end

        % Adaptive lambda from robust quantile of kappa^2
        k2 = cond_vals(finite_mask).^2;
        try
            qv = quantile(k2, opts.q);
        catch
            qv = prctile(k2, opts.q*100);
        end
        if ~isfinite(qv) || qv <= 0
            medk2 = median(k2(~isnan(k2) & isfinite(k2)));
            if ~isfinite(medk2) || medk2 <= 0, medk2 = 1; end
            qv = medk2;
        end
        lambda_g = opts.alpha * (d / qv);
        lambda_g = min(max(lambda_g, opts.lambda_bounds(1)), opts.lambda_bounds(2));

        % Final scores; pick argmax (no cross-sector constraints)
        scores = pres_vals - lambda_g * (cond_vals.^2);
        [best_score, wbest] = max(scores);
        selected_pairs{g} = sort([wrap(wbest), wrap(wbest+1)]).';

        meta.lambda(g)     = lambda_g;
        meta.best_score(g) = best_score;
        meta.best_kappa(g) = cond_vals(wbest);
    end
end