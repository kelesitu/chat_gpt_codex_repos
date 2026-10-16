function [minDB, wMin] = svMinimumBefore(w, Tgrid, Tfun, wLast)
% Her tekil deger icin w(1) < frekans < wLast(k) araligindaki
% en derin YEREL cukuru bulur. Cukur yoksa iki cikti da NaN olur.
%
% w      : Artan pozitif frekanslar [rad/s]
% Tgrid  : Kompleks frekans cevaplari [ny x nu x numel(w)]
% Tfun   : Tek frekansta T matrisini veren fonksiyon
% wLast  : Her tekil degerin ust frekans siniri [rad/s]
% minDB  : Secilen cukurun genligi [dB]
% wMin   : Secilen cukurun frekansi [rad/s]
%
% Gridde gorunmeyen dar cukurlar ve duz minimum bolgeleri kacabilir.
% Tgrid ve Tfun ayni sistemi temsil etmelidir.

    validateattributes(w, {'numeric'}, ...
        {'real', 'finite', 'positive', 'vector', 'nonempty'});
    validateattributes(wLast, {'numeric'}, {'real', 'vector'});
    assert(isa(Tfun, 'function_handle'), 'Tfun fonksiyon olmali.');

    w = w(:).';
    wLast = wLast(:);
    nSV = min(size(Tgrid, 1), size(Tgrid, 2));
    assert(all(diff(w) > 0), 'w kesin artan olmali.');
    assert(size(Tgrid, 3) == numel(w), 'Tgrid ve w boyutlari uyusmuyor.');
    assert(all(isfinite(Tgrid(:))), 'Tgrid sonlu degerler icermeli.');
    assert(numel(wLast) == nSV, 'Her tekil deger icin bir sinir gerekli.');

    % 1. Griddeki tekil degerleri bir kez hesapla.
    sigma = zeros(nSV, numel(w));
    for j = 1:numel(w)
        sigma(:, j) = svd(Tgrid(:, :, j));
    end

    minDB = nan(nSV, 1);
    wMin = nan(nSV, 1);
    options = optimset('Display', 'off', 'TolX', 1e-10);

    for k = 1:nSV
        limit = wLast(k);
        if ~isfinite(limit)
            continue
        end
        assert(limit >= w(1) && limit <= w(end), ...
            'Frekans siniri grid araliginda olmali.');

        % 2. Araligi sinir frekansinda kes; yerel cukurlari belirle.
        response = @(x) singularValue(Tfun, exp(x), k);
        inside = w < limit;
        x = log([w(inside), limit]);
        y = [sigma(k, inside), response(log(limit))];

        slope = diff(y);
        dips = find(slope(1:end-1) <= 0 & slope(2:end) >= 0 & ...
                   (slope(1:end-1) < 0 | slope(2:end) > 0)) + 1;

        % 3. Cukurlari hassaslastir ve en derin olanini sec.
        best = Inf;
        for j = dips
            [value, location] = refineDip(response, ...
                x(j-1:j+1), y(j-1:j+1), options);
            frequency = exp(location);

            if value < best && frequency > w(1) && frequency < limit
                best = value;
                minDB(k) = 20*log10(value);
                wMin(k) = frequency;
            end
        end
    end
end

function [value, location] = refineDip(response, x, y, options)
% Sinira kayan veya iki yaninda yukselis olmayan sonucu ele.
    [location, value, flag] = fminbnd(response, x(1), x(3), options);
    if flag <= 0 || location <= x(1) || location >= x(3) || ...
            value >= y(1) || value >= y(3)
        value = Inf;
        return
    end

    % Gercek yerel cukur kontrolu: iki yakin noktada genlik daha buyuk olmali.
    step = min(1e-3*(x(3)-x(1)), ...
               0.25*min(location-x(1), x(3)-location));
    if step <= 0 || location-step == location || location+step == location
        value = Inf;
        return
    end
    neighbors = [response(location-step), response(location+step)];
    tolerance = 32*eps(max(abs([value, neighbors])));
    if ~all(neighbors-value > tolerance)
        value = Inf;
    end
end

function value = singularValue(Tfun, frequency, k)
    values = svd(Tfun(frequency));
    value = values(k);
end
