function beam_groups = sectorizeDftBeams(mu, M, W, zeta)
% Sectorize DFT beamspace using sliding-window sectorization method
%
% Inputs:
%   mu    - spatial frequencies in radians (between -pi and pi), size [L x 1]
%   M     - number of antennas (equals number of DFT beams)
%   W     - window size (number of DFT beams per sector)
%   zeta  - max beam-index offset allowed from window center (symmetry rule)
%
% Output:
%   beam_groups - cell array containing beam index vectors for each SoI

    % Normalize mu to [0, 2*pi)
    mu = mod(mu, 2*pi);
    
    % Sort in ascending order
    mu = sort(mu(:));  % ensure column vector

    % Discretized beam centers
    L = numel(mu);
    unclassified = true(L,1);
    beam_groups = {};

    while any(unclassified)
        % Candidate SoIs starting from all possible beam indices
        best_group = [];
        best_score = -Inf;

        for k0 = 1:M
            % Define window [k0, ..., k0+W-1] modulo M
            beams = mod((k0-1):(k0+W-2), M) + 1;

            % Window center
            mu_center = mod((2*pi/M) * mod(k0 + W/2 - 1, M), 2*pi);

            % Compute distance from center (handle wrap-around)
            dists = abs(angle(exp(1j*(mu(unclassified) - mu_center))));
            in_window = dists <= (zeta * pi / M);

            if sum(in_window) > best_score
                best_score = sum(in_window);
                best_group = beams;
                best_unclassified = find(unclassified);
                group_members = best_unclassified(in_window);
            end
        end

        if isempty(best_group)
            warning('Could not assign all frequencies. Increase zeta.');
            break;
        end

        beam_groups{end+1} = best_group; %#ok<AGROW>
        unclassified(group_members) = false;
    end
end