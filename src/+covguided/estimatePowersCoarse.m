function [p_est, sigma2_est, mu_sorted, R_est] = estimatePowersCoarse(R_hat, maskIdx, mu_hat, M)
    % Inputs:
    %   R_hat  : [NRF x NRF] covariance matrix
    %   mu_hat : [d x 1] estimated spatial frequencies
    %
    % Outputs:
    %   p_est      : [d x 1] estimated powers
    %   sigma2_est : scalar noise power
    
    %%
    d = length(mu_hat);
    NRF = size(R_hat, 1);

    % Build A on full aperture, then index rows:
    A = @(mu) exp(1i * (0:M-1).' * mu.'); 
    Afull  = A(mu_hat);                   % [M x d]
    A_mask = Afull(maskIdx, :);           % [NRF x d]
    
    % --- Quadratic term Q and linear term b without forming B
    G    = A_mask' * A_mask;
    Q11  = abs(G).^2;
    q    = real(sum(abs(A_mask).^2,1)).';
    Q    = [Q11, q; q.', NRF];
    Q    = real(0.5*(Q+Q.'));         % exact symmetric & real
    
    bsig = real(sum(conj(A_mask).*(R_hat*A_mask),1)).';
    b    = [bsig; real(trace(R_hat))];
    
    % Tiny ridge if needed
    epsQ = 1e-12 * trace(Q)/size(Q,1);
    Q    = Q + epsQ * eye(d+1);

    % --- Solve QP (x >= 0)
    opts = optimoptions('quadprog','Display','off','Algorithm','interior-point-convex');
    x_est = quadprog(2*Q, -2*b, [], [], [], [], zeros(d+1,1), [], [], opts);
    
    p_est     = x_est(1:d);
    sigma2_est = x_est(end);


    % Step 7: Reconstruct signal covariance
    [p_sorted, idx_p] = sort(p_est,'descend');
    mu_sorted = mu_hat(idx_p);
    R = Afull(:, idx_p) * diag(p_sorted) * Afull(:, idx_p)';

    %% Toeplitz plus PSD Projection
    R_est = covguided.toeplitzProjection(R);

end
