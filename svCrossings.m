function [wFirst, wLast, wAll, wLastMinus3] = ...
    svCrossings(w, Hgrid, Hfun, levelDB)
%SVCROSSINGS Crossings of sorted singular values.
% First three outputs refer to levelDB (default 0 dB).
% Fourth output contains the last -3 dB crossing of each singular value.
% w: positive increasing frequencies [rad/s].
% Hgrid: ny-by-nu-by-numel(w) complex responses.
% Hfun(wi): response matrix at one frequency, consistent with Hgrid.
% Missing crossings return NaN. wAll{k} contains all detected crossings.
% Hidden crossing pairs and tangencies may be missed by the initial grid.

    if nargin < 4
        levelDB = 0;
    end
    validateattributes(levelDB, {'numeric'}, ...
        {'real', 'finite', 'scalar'});
    validateattributes(w, {'numeric'}, ...
        {'real', 'finite', 'positive', 'vector', 'nonempty'});
    assert(isa(Hfun, 'function_handle'), 'Hfun must be a function handle.');

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
    wLast = nan(nSV, 1);
    wAll = cell(nSV, 1);
    wLastMinus3 = nan(nSV, 1);

    for k = 1:nSV
        wc = findCrossings(x, w, S(k,:), Hfun, k, levelDB, opts);
        wAll{k} = wc(:);
        if ~isempty(wc)
            wFirst(k) = wc(1);
            wLast(k) = wc(end);
        end

        if nargout >= 4
            if levelDB == -3
                wc3 = wc;
            else
                wc3 = findCrossings(x, w, S(k,:), Hfun, k, -3, opts);
            end
            if ~isempty(wc3)
                wLastMinus3(k) = wc3(end);
            end
        end
    end
end

function wc = findCrossings(x, w, sGrid, Hfun, k, levelDB, opts)
    target = 10^(levelDB/20);
    assert(isfinite(target) && target > 0, ...
        'Target level is outside the supported numerical range.');
    g = sGrid - target;
    wc = w(g == 0);
    idx = find( ...
        (g(1:end-1) < 0 & g(2:end) > 0) | ...
        (g(1:end-1) > 0 & g(2:end) < 0));

    f = @(z) svResidual(z, k, Hfun, target);
    for j = idx
        xr = fzero(f, [x(j), x(j+1)], opts);
        % Preserve original interval bounds against exp/log roundoff.
        wc(end+1) = min(w(j+1), max(w(j), exp(xr))); %#ok<AGROW>
    end
    wc = unique(wc);
end

function y = svResidual(x, k, Hfun, target)
    s = svd(Hfun(exp(x)));
    y = s(k) - target;
end
