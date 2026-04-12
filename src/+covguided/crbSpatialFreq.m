% src/+covguided/crbSpatialFreq.m
function c = crbSpatialFreq(M, N, d, asnrDb)
    % Stochastic CRB bound approximation for ULA, d equal-power sources,
    % returns sqrt(CRB) per ASNR point (rad).
    snr = 10.^(asnrDb/10);
    c = sqrt( 6 ./ (snr * N * M * (M^2-1)) ) * sqrt(d);
end