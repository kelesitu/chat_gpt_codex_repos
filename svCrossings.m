function [wFirst, wLast, wAll, wLastMinus3] = ...
    svCrossings(w, Hgrid, Hfun, levelDB)
%SVCROSSINGS Tekil degerlerin seviye kesisimlerini bulur.
% Ilk uc cikti levelDB (varsayilan 0 dB), son cikti -3 dB icindir.
% w: artan pozitif frekanslar [rad/s].
% Hgrid: ny x nu x numel(w) kompleks frekans cevaplari.
% Hfun(wi): tek frekansta ayni sistemin cevap matrisi.
% Kesisim yoksa NaN; wAll{k} tum TESPIT EDILEN kesisimlerdir.
%
% Ilk kesisimden onceki cukurlar ve son kesisimden sonraki tepeler,
% en yakin adaydan baslayarak hassaslastirilir. Gridde ilk basta
% kesisim yoksa butun cukur ve tepe adaylari incelenir.
% Gridde hic izi olmayan dar salinimlar yine kacabilir.

    if nargin < 4
        levelDB = 0;
    end
    validateattributes(levelDB, {'numeric'}, ...
        {'real', 'finite', 'scalar'});
    validateattributes(w, {'numeric'}, ...
        {'real', 'finite', 'positive', 'vector', 'nonempty'});
    assert(isa(Hfun, 'function_handle'), 'Hfun fonksiyon olmali.');

    w = w(:).';
    assert(all(diff(w) > 0), 'w kesin artan olmali.');
    assert(size(Hgrid, 3) == numel(w), 'Hgrid ve w boyutlari uyusmuyor.');
    assert(all(isfinite(Hgrid(:))), 'Hgrid sonlu degerler icermeli.');

    nSV = min(size(Hgrid, 1), size(Hgrid, 2));
    sigma = zeros(nSV, numel(w));
    for j = 1:numel(w)
        sigma(:, j) = svd(Hgrid(:, :, j));
    end

    x = log(w);
    options = optimset('Display', 'off', 'TolX', 1e-10);
    wFirst = nan(nSV, 1);
    wLast = nan(nSV, 1);
    wAll = cell(nSV, 1);
    wLastMinus3 = nan(nSV, 1);

    for k = 1:nSV
        wc = findCrossings(x, w, sigma(k,:), Hfun, k, levelDB, options);
        wAll{k} = wc(:);
        if ~isempty(wc)
            wFirst(k) = wc(1);
            wLast(k) = wc(end);
        end

        if nargout >= 4
            if levelDB == -3
                wc3 = wc;
            else
                wc3 = findCrossings(x, w, sigma(k,:), Hfun, k, -3, options);
            end
            if ~isempty(wc3)
                wLastMinus3(k) = wc3(end);
            end
        end
    end
end

function wc = findCrossings(x, w, sigma, Hfun, k, levelDB, options)
    target = 10^(levelDB/20);
    assert(isfinite(target) && target > 0, ...
        'Hedef seviye sayisal araligin disinda.');
    g = sigma - target;
    response = @(z) svResidual(z, k, Hfun, target);

    % 1. Gridde gorulen kesisimleri hassaslastir.
    roots = x(g == 0);
    intervals = find((g(1:end-1) < 0 & g(2:end) > 0) | ...
                     (g(1:end-1) > 0 & g(2:end) < 0));
    for j = intervals
        roots = [roots, bracketRoot(response, x(j), x(j+1), ...
                                   g(j), g(j+1), options)]; %#ok<AGROW>
    end
    roots = sort(roots);

    % 2. Egim degisimlerinden aday indekslerini kaydet.
    slope = diff(sigma);
    minima = find(slope(1:end-1) <= 0 & slope(2:end) >= 0 & ...
                 (slope(1:end-1) < 0 | slope(2:end) > 0)) + 1;
    maxima = find(slope(1:end-1) >= 0 & slope(2:end) <= 0 & ...
                 (slope(1:end-1) > 0 | slope(2:end) < 0)) + 1;

    % Sinirlari ilk tespit edilen kesisimlere gore belirle.
    if ~isempty(roots)
        minima = minima(x(minima) < roots(1));
        maxima = maxima(x(maxima) > roots(end));
    end

    % 3. Soldaki cukurlari en yakin olandan baslayarak tara.
    for j = fliplr(minima)
        extra = refineCandidate(response, x, g, j, +1, options);
        roots = sort([roots, extra]); %#ok<AGROW>
        % roots(1), yeni ilk kesisimdir. Daha soldaki adaylara devam et.
    end

    % 4. Sagdaki tepeleri en yakin olandan baslayarak tara.
    for j = maxima
        extra = refineCandidate(response, x, g, j, -1, options);
        roots = sort([roots, extra]); %#ok<AGROW>
        % roots(end), yeni son kesisimdir. Daha sagdaki adaylara devam et.
    end

    % Sayisal olarak ayni kokleri birlestir.
    wc = [];
    if ~isempty(roots)
        tolerance = max(options.TolX, 32*eps(max(1, max(abs(roots)))));
        roots = roots([true, diff(roots) > tolerance]);
        wc = min(w(end), max(w(1), exp(roots)));
    end
end

function roots = refineCandidate(response, x, g, j, direction, options)
% direction=+1: minimum, direction=-1: maksimum ara.
    [extremum, value, flag] = fminbnd( ...
        @(z) direction*response(z), x(j-1), x(j+1), options);
    assert(flag > 0, 'Aday hassaslastirmasi yakinsamadi.');
    value = direction*value;

    % Ekstremum seviyenin diger tarafindaysa iki yeni kok bulunabilir.
    left = bracketRoot(response, x(j-1), extremum, ...
                       g(j-1), value, options);
    right = bracketRoot(response, extremum, x(j+1), ...
                        value, g(j+1), options);
    roots = [left, right];
end

function root = bracketRoot(response, a, b, ga, gb, options)
% Yalnizca sifir deger veya isaret degisimi varsa kok ara.
    root = [];
    if ga == 0
        root(end+1) = a;
    end
    if gb == 0
        root(end+1) = b;
    end
    if (ga < 0 && gb > 0) || (ga > 0 && gb < 0)
        [value, ~, flag] = fzero(response, [a, b], options);
        assert(flag > 0, 'Kesisim aramasi yakinsamadi.');
        root(end+1) = value;
    end
end

function value = svResidual(x, k, Hfun, target)
    values = svd(Hfun(exp(x)));
    value = values(k) - target;
end
