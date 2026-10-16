function [minDB, wMin] = svMinimumBefore(w, Tgrid, Tfun, wLast)
%SVMINIMUMBEFORE Deepest interior local minimum before each frequency limit.
% w: positive increasing frequencies [rad/s].
% Tgrid: ny-by-nu-by-numel(w) complex responses.
% Tfun(wi): response matrix at one frequency, consistent with Tgrid.
% wLast(k): strict upper frequency limit for singular value k.
% Searches w(1) < frequency < wLast(k), excluding both endpoints.
% Only numerically resolved, strict local minima are accepted.
% If none are found, or the limit is nonfinite, both outputs are NaN.
% Candidate minima visible on the grid are refined in log frequency.
% Hidden narrow minima and flat minimum intervals may be missed.
% A zero minimum returns -Inf dB.
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

        % Endpoints bracket candidates but are never minimum candidates.
        idx = find( ...
            sk(2:end-1) <= sk(1:end-2) & ...
            sk(2:end-1) <= sk(3:end) & ...
           (sk(2:end-1) < sk(1:end-2) | ...
            sk(2:end-1) < sk(3:end))) + 1;

        bestS = Inf;
        bestW = NaN;
        f = @(x) svAt(x, k, Tfun);

        for j = idx
            a = xk(j-1);
            b = xk(j+1);
            [xm, sm, exitflag] = fminbnd(f, a, b, opts);

            % Reject a solver result that lies at a bracket boundary or
            % fails to form a valley relative to both bracket endpoints.
            if exitflag <= 0 || ~(a < xm && xm < b) || ...
                    ~(sm < sk(j-1) && sm < sk(j+1))
                continue
            end

            % Verify an actual local dip on both sides. No fixed frequency
            % exclusion band: genuine dips close to wLast remain eligible.
            h = min(1e-3*(b-a), 0.25*min(xm-a, b-xm));
            if h <= 0 || xm-h == xm || xm+h == xm
                continue
            end
            left = f(xm-h);
            right = f(xm+h);
            tolS = 32*eps(max([abs(sm), abs(left), abs(right)]));
            if ~(left-sm > tolS && right-sm > tolS)
                continue
            end

            wm = exp(xm);
            if w(1) < wm && wm < upper && sm < bestS
                bestS = sm;
                bestW = wm;
            end
        end

        if isfinite(bestW)
            minDB(k) = 20*log10(bestS);
            wMin(k) = bestW;
        end
    end
end

function value = svAt(x, k, Tfun)
    s = svd(Tfun(exp(x)));
    value = s(k);
end
