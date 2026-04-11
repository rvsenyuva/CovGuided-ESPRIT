% ---------- FUNCTION 1: Subspace Angle Errors (Haardt + Principal Angle) ----------
function theta_haardt = computeSubspaceAngleMetrics(U_true, U_est, d)
% Inputs:
%   U_true: MxM true covariance matrix
%   U_est : MxM estimated covariance matrix
%   d     : number of signal subspace dimensions
% Outputs:
%   theta_haardt   : Haardt-style subspace angle error (in degrees)
%   theta_principal: Principal angle-based RMS error (in degrees)

% Extract d-dimensional signal subspaces
[U1, ~] = eigs(U_true, d);
[U2, ~] = eigs(U_est, d);

% Haardt-style metric
P1 = U1 * U1';
P2 = U2 * U2';
theta_haardt = acosd( real( trace(P1 * P2) / d ) );

% Principal angle-based RMS
%s = svd(U1' * U2);
%theta_principal = sqrt(mean((acosd(min(s,1))).^2));
end