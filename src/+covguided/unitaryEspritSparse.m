function mu_UE = unitaryEspritSparse(Y, beamCols, d, M)
% Unitary ESPRIT using selected beamCols from phase-shifted DFT basis
% Inputs:
%   Y        : [NRF x N] beamspace data matrix
%   beamCols : [L x 1] 1-based indices of selected beamforming vectors
%   d        : Model order (number of sources)
%   M        : Total number of antennas (for canonical basis generation)
% Output:
%   mu_UE    : [d_eff x 1] Estimated spatial frequencies (radians)

    % === Step 0: Real FBA augmentation for Unitary ESPRIT
    Y_UE = sqrt(2) * [real(Y),  imag(Y)];

    % === Step 1: Real-valued SVD subspace extraction
    [Us_UE, ~, ~] = svds(Y_UE, d);                  % [NRF x <=d]
    d_eff         = size(Us_UE, 2);                 % actual #cols we have

    % === Step 2: Canonical cosine/sine transform matrices
    cosvec = cos((0:M-1) * pi/M);  % 1×M
    sinvec = sin((0:M-1) * pi/M);  % 1×M
   
    % === Step 3: Selection rows that form a valid forward shift
    cols      = unique(beamCols(:));
    next_idx  = mod(cols, M) + 1;
    beamRows  = cols(ismember(next_idx, cols));
    beamNext  = mod(beamRows, M) + 1;  % i → i+1 mod M
    L         = numel(beamRows);
    
    % === Step 4: Populate E1, E2 using projected subspace rows

    E1 = zeros(L, d_eff);
    E2 = zeros(L, d_eff);

    for l = 1:L
        i   = beamRows(l);
        ip1 = beamNext(l);

        idx_i   = find(beamCols == i,   1);
        idx_ip1 = find(beamCols == ip1, 1);

        % (Paranoid) skip if either is missing
        %if isempty(idx_i) || isempty(idx_ip1), continue; end

        % Canonical weights (J1, J2) — wraparound handling
        if i == M && ip1 == 1
            c1 = (-1)^M; s1 = 0;  % Special wrap case
        else
            c1 = cosvec(ip1);
            s1 = sinvec(ip1);
        end
        c0 = cosvec(i);
        s0 = sinvec(i);

        % Linear projection of subspace rows (1×d_eff each)
        E1(l, :) = c0 * Us_UE(idx_i, :) + c1 * Us_UE(idx_ip1, :);
        E2(l, :) = s0 * Us_UE(idx_i, :) + s1 * Us_UE(idx_ip1, :);
    end

    % Drop any zero rows (if a pair was skipped)
    nz = any(E1,2) | any(E2,2);
    E1 = E1(nz,:);  E2 = E2(nz,:);
    L  = size(E1,1);

    sE  = svd(E1,'econ');
    condE = sE(1)/max(sE(end), eps);
    ill = (L < d_eff) || (condE > 1e6);   % threshold: tune 1e5–1e8

    % === Step 5: Solve the real-valued shift-invariance equation
    if ill
        % ---- UNDERDETERMINED: use ridge regularization ----
        % Auto lambda scaled to data (tune 1e-3..1e-1 if needed)
        sE  = svd(E1,'econ');
        lam = 1e-2 * (median(sE)^2 + realmin);

        S1 = E1.'*E1;
        S2 = E1.'*E2;

        % SPD solve via Cholesky
        R   = chol(S1 + lam*eye(d_eff), 'lower');
        Phi = R'\(R\S2);
    else
        % ---- DETERMINED/OVERDETERMINED: your original path ----
        % Phi = lsqminnorm(E1, E2);
        % well-conditioned normal equations
        S1 = E1.'*E1;  
        S2 = E1.'*E2;
        
        % SPD solve via Cholesky
        R  = chol(S1 + 1e-12*trace(S1)/d_eff*eye(d_eff), 'lower');
        Phi = R'\(R\S2);
    end

    % === Step 6: Frequency extraction (keep your original style)
    [~, T] = schur(Phi, 'real');
    mu_UE = 2 * atan(diag(T));  % Unitary ESPRIT frequency estimates
end