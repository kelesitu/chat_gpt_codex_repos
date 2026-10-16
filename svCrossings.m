function [wFirst, wLast, wAll] = svCrossings(w, Hgrid, Hfun, levelDB)
%SVCROSSINGS Find first/last level crossings of each sorted singular value.
% w: positive increasing frequencies [rad/s].
% Hgrid: ny-by-nu-by-numel(w) complex frequency response.
% Hfun(wi): complex response matrix at one frequency.
% levelDB: target magnitude in dB (default: 0).
% Outputs: first/last crossings [rad/s]; NaN when none are detected.
% wAll{k} contains all detected crossings of the kth largest singular value.
% Only sign changes and exact grid hits are detected. Hidden pairs of
% crossings and tangencies can be missed by the initial grid.

    if nargin < 4
        levelDB = 0;
    end
    validateattributes(levelDB, {'numeric'}, ...
        {'real', 'finite', 'scalar'});
    validateattributes(w, {'numeric'}, ...
        {'real', 'finite', 'positive', 'vector', 'nonempty'});
    assert(isa(Hfun, 'function_handle'), 'Hfun must be a function handle.');

    target = 10^(levelDB/20);
    assert(isfinite(target) && target > 0, ...
        'levelDB is outside the supported numerical range.');
    w = w(:).';
    assert(all(diff(w) > 0), 'w must be strictly increasing.');
    assert(size(Hgrid, 3) == numel(w), ...
        'Hgrid and w sizes must agree.');
    assert(all(isfinite(Hgrid(:))), 'Hgrid must contain finite values.');

    nSV = min(size(Hgrid, 1), size(Hgrid, 2));
    S = zeros(nSV, numel(w));
    for j = 1:numel(w)
        S(:, j) = svd(Hgrid(:, :, j));
    end

    x = log(w);
    opts = optimset('Display', 'off', 'TolX', 1e-10);
    wFirst = nan(nSV, 1);
    wLast  = nan(nSV, 1);
    wAll   = cell(nSV, 1);

    for k = 1:nSV
        g = S(k, :) - target;
        rootsX = x(g == 0);
        idx = find( ...
            (g(1:end-1) < 0 & g(2:end) > 0) | ...
            (g(1:end-1) > 0 & g(2:end) < 0));

        f = @(z) svResidual(z, k, Hfun, target);
        for j = idx
            rootsX(end+1) = fzero(f, [x(j), x(j+1)], opts); %#ok<AGROW>
        end

        wc = exp(unique(rootsX));
        wAll{k} = wc(:);
        if ~isempty(wc)
            wFirst(k) = wc(1);
            wLast(k)  = wc(end);
        end
    end
end

function y = svResidual(x, k, Hfun, target)
    s = svd(Hfun(exp(x)));
    y = s(k) - target;
end
