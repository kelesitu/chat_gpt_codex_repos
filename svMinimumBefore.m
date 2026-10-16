function [minDB, wMin] = svMinimumBefore(w, Tgrid, Tfun, wLast)
%SVMINIMUMBEFORE Minimum sorted singular values up to given frequencies.
% w: positive increasing frequencies [rad/s].
% Tgrid: ny-by-nu-by-numel(w) complex responses.
% Tfun(wi): response matrix at one frequency, consistent with Tgrid.
% wLast(k): inclusive upper frequency limit for singular value k.
% Search starts at w(1). Nonfinite limits return NaN.
% Local minima visible on the grid are refined in log frequency.
% Narrow hidden minima may be missed. A zero minimum returns -Inf dB.
%
% Example:
% [~,~,~,wLastMinus3] = svCrossings(w,Lgrid,Lfun,0);
% [minDB,wMin] = svMinimumBefore(w,Tgrid,Tfun,wLastMinus3);

    validateattributes(w, {'numeric'}, ...
        {'real', 'finite', 'positive', 'vector', 'nonempty'});
    validateattributes(wLast, {'numeric'}, {'real', 'vector'});
    assert(isa(Tfun, 'function_handle'), 'Tfun must be a function handle.');
    w = w(:).';
    wLast = wLast(:);

    assert(all(diff(w) > 0), 'w must be strictly increasing.');
    assert(size(Tgrid, 3) == numel(w), ...
        'Tgrid and w sizes must agree.');
    assert(all(isfinite(Tgrid(:))), 'Tgrid must contain finite values.');

    nSV = min(size(Tgrid, 1), size(Tgrid, 2));
    assert(numel(wLast) == nSV, ...
        'Provide one upper frequency limit per singular value.');

    S = zeros(nSV, numel(w));
    for j = 1:numel(w)
        S(:, j) = svd(Tgrid(:, :, j));
    end

    minDB = nan(nSV, 1);
    wMin = nan(nSV, 1);
    opts = optimset('Display', 'off', 'TolX', 1e-10);

    for k = 1:nSV
        upper = wLast(k);
        if ~isfinite(upper)
            continue
        end
        assert(upper >= w(1) && upper <= w(end), ...
            'Upper frequency limit must lie within the grid range.');

        inside = w < upper;
        wk = [w(inside), upper];
        sk = [S(k, inside), svAt(log(upper), k, Tfun)];
        xk = log(wk);

        [bestS, jBest] = min(sk);
        bestW = wk(jBest);

        idx = find( ...
            sk(2:end-1) <= sk(1:end-2) & ...
            sk(2:end-1) <= sk(3:end) & ...
           (sk(2:end-1) < sk(1:end-2) | ...
            sk(2:end-1) < sk(3:end))) + 1;

        f = @(x) svAt(x, k, Tfun);
        for j = idx
            [xm, sm] = fminbnd(f, xk(j-1), xk(j+1), opts);
            if sm < bestS
                bestS = sm;
                bestW = min(wk(j+1), max(wk(j-1), exp(xm)));
            end
        end

        minDB(k) = 20*log10(bestS);
        wMin(k) = bestW;
    end
end

function value = svAt(x, k, Tfun)
    s = svd(Tfun(exp(x)));
    value = s(k);
end
