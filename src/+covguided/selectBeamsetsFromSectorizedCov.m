function [selected_sets, meta] = selectBeamsetsFromSectorizedCov(R_model, beam_groups, K, opts)
% Generalizes select_adjacent_pairs_from_sectorized_cov to K >= 2 beams/sector.
% Fast path: precomputes full DFT bank and uses capture = tr((G+γI)^{-1} (B^H Rs B)).
% Optional pruning: keep only windows whose center is among top-q beam energies in beamspace.
%
% Inputs:
%   R_model     : [M x M] Hermitian PSD covariance (e.g., from power fitting)
%   beam_groups : 1xG cell, each a vector of candidate beam indices (SoIs, 1..M)
%   K           : scalar or [G x 1] per-sector K_g, K_g >= 2
%   opts        : struct with fields
%                 .eig_tol_scale   (default 1e-6)
%                 .tikhonov_gamma  (default 0)
%                 .alpha           (default 0.5)  % conditioning penalty weight
%                 .normalize_score (default true)
%                 .prune_topq      (default 0)    % 0 disables pruning; try 3 or 4 for speed
%
% Outputs:
%   selected_sets : 1xG cell, each a 1xK_g vector (selected beam indices for sector g)
%   meta          : .best_score, .best_kappa, .Kg, .num_shifts, .d, .pruned_windows

    if nargin < 4, opts = struct(); end
    if ~isfield(opts,'eig_tol_scale'),   opts.eig_tol_scale   = 1e-6; end
    if ~isfield(opts,'tikhonov_gamma'),  opts.tikhonov_gamma  = 0;    end
    if ~isfield(opts,'alpha'),           opts.alpha           = 0.5;  end
    if ~isfield(opts,'normalize_score'), opts.normalize_score = true; end
    if ~isfield(opts,'prune_topq'),      opts.prune_topq      = 0;    end

    Gs = numel(beam_groups);
    M  = size(R_model,1);
    if isscalar(K), K = repmat(K, Gs, 1); end

    % --- Soft-thresholded signal covariance proxy Rs ---
    R = (R_model + R_model')/2;
    [U, Svals] = eig(R, 'vector'); [Svals, idx] = sort(Svals, 'descend'); U = U(:, idx);
    tol = opts.eig_tol_scale * norm(R, 'fro');
    d  = sum(Svals > max(tol,0));  if d < 1, error('Detected signal rank is zero.'); end
    tau = 0; if d < numel(Svals), tau = median(Svals(d+1:end)); end
    lam = max(Svals(1:d) - tau, 0);
    Rs  = U(:,1:d) * diag(lam) * U(:,1:d)';                % "denoised" signal cov
    if opts.normalize_score, Rs = Rs / max(real(trace(Rs)), eps); end

    % --- Precompute full DFT-like bank and helpers ---
    n = (0:M-1).';
    allbins0 = 0:M-1;
    Bfull  = exp(1i*2*pi*n*(allbins0/M)) .* exp(-1i*(M-1)*pi*(allbins0/M));   % [M x M]
    RsBfull = Rs * Bfull;                                                     % [M x M]
    % Beamspace diagonal energy for pruning: diag(B^H Rs B)
    Rb_diag = real(sum(conj(Bfull) .* RsBfull, 1));                           % 1 x M

    selected_sets = cell(Gs,1);
    meta = struct('best_score', nan(Gs,1), 'best_kappa', nan(Gs,1), ...
                  'Kg', nan(Gs,1), 'num_shifts', nan(Gs,1), ...
                  'd', d, 'pruned_windows', nan(Gs,1));

    for g = 1:Gs
        cand = beam_groups{g}(:).';               % 1..M
        Lg   = numel(cand);
        Kg   = K(g);
        if Kg < 2, error('K must be >= 2.'); end
        if Kg > Lg
            warning('K=%d exceeds sector size=%d; clipping.', Kg, Lg);
            Kg = Lg;
        end

        % All possible starts of a length-Kg contiguous window
        starts = 1:(Lg - Kg + 1);

        % --- Prune by top-q center beams (optional) ---
        if opts.prune_topq > 0 && Lg > Kg
            q = min(opts.prune_topq, Lg);
            % top-q energies within this sector
            [~, ordE] = sort(Rb_diag(cand), 'descend');
            topCenters = sort(ordE(1:q));   % positions within 1..Lg
            half = (Kg-1)/2;
            keep = false(size(starts));
            for cpos = topCenters
                % Include all windows whose center overlaps this center position.
                smin = max(1, ceil(cpos - half));
                smax = min(Lg - Kg + 1, floor(cpos + half));
                keep(smin:smax) = true;
            end
            starts = starts(keep);
        end
        meta.pruned_windows(g) = (Lg - Kg + 1) - numel(starts);

        bestScore = -Inf; bestSet = []; bestKappa2 = Inf; bestCenter = Inf;

        for s = starts
            beams = cand(s:(s+Kg-1));                 % actual beam indices (1..M)
            Bsel  = Bfull(:, beams);                  % [M x Kg]
            % Small KxK Gram and Rs-projected Gram
            Gg = Bsel' * Bsel;                        % K x K
            Sv = svd(Gg);
            gamma = max(opts.tikhonov_gamma, 1e-6 * mean(Sv));
            invG = (Gg + gamma * eye(Kg)) \ eye(Kg);  % (Gg + γI)^{-1}

            S = Bsel' * RsBfull(:, beams);           % K x K equals B^H Rs B
            capture = real(trace(invG * S));         % tr((G+γI)^{-1} (B^H Rs B))
            kappa2  = cond(Gg)^2;
            score   = capture / (1 + opts.alpha * kappa2);

            % Tie-break: lower kappa, then more central window in the sector
            centerPos = s + (Kg-1)/2;                % fractional center in [1..Lg]
            choose = (score > bestScore) ...
                  || (abs(score - bestScore) <= 1e-12 && (kappa2 < bestKappa2 ...
                  || (abs(kappa2 - bestKappa2) <= 1e-12 && centerPos < bestCenter)));
            if choose
                bestScore  = score;
                bestSet    = beams;
                bestKappa2 = kappa2;
                bestCenter = centerPos;
            end
        end

        selected_sets{g}    = bestSet;
        meta.best_score(g)  = bestScore;
        meta.best_kappa(g)  = sqrt(bestKappa2);
        meta.Kg(g)          = Kg;
        meta.num_shifts(g)  = Kg - 1;
    end
end