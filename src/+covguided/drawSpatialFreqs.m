% src/+covguided/drawSpatialFreqs.m
function mu = drawSpatialFreqs(d, M, s)
    minSep = 2*pi/M;                        % one DFT bin
    while true
        mu = -pi + 2*pi*rand(s, d, 1);
        if d == 1 || min(diff(sort(mu))) >= minSep, return; end
    end
end