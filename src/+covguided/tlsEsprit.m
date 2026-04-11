function mu_est = tlsEsprit(R, d)
% TLS-ESPRIT with basic safeguards.
% R     : NRF x NRF covariance (FBA-averaged)
% J1,J2 : selection matrices for given shift Delta
% d     : model order
% Delta : shift (default 1)
%
% mu_est: d x 1 spatial frequencies in radians

NRF = size(R, 1);

% 2) Signal subspace (top-d eigenvectors)
% eigs is fine for larger NRF, eig is fine for small NRF; keep eig for determinism
[U, S] = eig(R, 'vector');
[~, idx] = sort(real(S), 'descend');
Us = U(:, idx(1:d));                      % NRF x d

% 3) Build TLS data matrix
E1 = Us(1:NRF-1, :);                       % ≡ J1*Us, with J1 = [I_{L-1}  0]
E2 = Us(2:NRF,   :);                       % ≡ J2*Us, with J2 = [0  I_{L-1}]
Z  = [E1, E2];                             % (NRF-Delta) x (2d)

% 4) TLS via SVD of Z
[~, ~, V] = svd(Z, 'econ');                % V: (2d x 2d)
V12 = V(1:d,      d+1:end);                % d x d
V22 = V(d+1:end,  d+1:end);                % d x d

% Condition check; regularize if necessary
rcondV22 = rcond(V22);
if rcondV22 < 1e-10
    % mild Tikhonov for numerical safety
    V22 = V22 + (1e-10 * norm(V22,'fro')/d) * eye(d);
end

% TLS solution: Phi = -V12 / V22
Phi = - V12 / V22;

% 5) Extract frequencies: eigenvalues of Phi ≈ exp(j*Delta*mu)
mu_est = angle(eig(Phi));

% Optional: sort by angle or by magnitude of eigenvalues closeness to unit circle
% [~, sidx] = sort(mu_est); mu_est = mu_est(sidx);

end