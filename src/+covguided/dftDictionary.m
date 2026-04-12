function B = dftDictionary(M)
%DFTDICTIONARY  Centered unitary DFT dictionary matching
%   phase_k .* fft(Y,M,1) in runMainExperiment.m.
%   Columns are DFT steering vectors post-multiplied by unit-modulus
%   phases so that B' * Y == phase_k .* fft(Y, M, 1).
    F       = fft(eye(M), M, 1);             % standard DFT matrix
    phase_k = exp(1i*(M-1)*pi*(0:M-1).'/M);  % centering phases
    B       = F' .* conj(phase_k).';         % M x M, columns unit-modulus
end