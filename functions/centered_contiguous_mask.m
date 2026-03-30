function idx = centered_contiguous_mask(M, NRF)
    % Preconditions for exact symmetry:
    % - if mod(M,2)==0, require mod(NRF,2)==0
    if mod(M,2)==0 && mod(NRF,2)==1
        error('For even M, NRF must be even to get contiguous centro-symmetry.');
    end
    i0  = floor((M - NRF)/2) + 1;
    idx = i0:(i0 + NRF - 1);
    
    
    W   = eye(M); 
    W   = W(idx,:);  % [NRF x M]

    % Sanity check (should be ~0):
    J_M   = flip(eye(M));
    J_NRF = flip(eye(NRF));
    assert(norm(J_NRF*W - W*J_M, 'fro') < 1e-12);
end
