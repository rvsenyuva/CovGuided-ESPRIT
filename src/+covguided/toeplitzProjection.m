function T_psd = toeplitzProjection(R)
% Projects a Hermitian matrix R onto the nearest Hermitian Toeplitz matrix
% Input:
%   R : [r x r] Hermitian sample covariance matrix
% Output:
%   T : [r x r] Hermitian Toeplitz matrix

r = size(R, 1);
T = zeros(r);

for k = -(r - 1):(r - 1)
    diag_vals = diag(R, k);
    avg_val = mean(diag_vals);
    T = T + diag(avg_val * ones(r - abs(k), 1), k);
end

% Enforce Hermitian symmetry explicitly
T_toeplitz = (T + T') / 2;

% PSD projection, i.e. eigenvalue clipping
[U, D] = eig(T_toeplitz);
d = real(diag(D));
d(d < 0) = 0;
T_psd = U * diag(d) * U';

end