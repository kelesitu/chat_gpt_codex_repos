function test_svCrossings
% Hazir transfer fonksiyonu veya Control System Toolbox gerektirmez.
% MATLAB'da bu dosyanin oldugu klasorde test_svCrossings calistirin.
    w = exp(0:12);
    width = 0.32;
    offset = width*sqrt(log(2));
    expected = [1.4-offset, 1.4+offset, 3.4-offset, 3.4+offset, 6, ...
                8.4-offset, 8.4+offset, 10.4-offset, 10.4+offset];

    % Gridde yalniz bir kesisim var; dort gizli cift bulunmali.
    Hfun = @(wi) exp(profile(log(wi)));
    Hgrid = sampleResponse(w, Hfun);
    g = reshape(Hgrid, 1, []) - 1;
    assert(sum(g == 0) + sum(g(1:end-1).*g(2:end) < 0) == 1);
    [first, last, all, last3] = svCrossings(w, Hgrid, Hfun, 0);
    assertRoots(all{1}, expected);
    assert(abs(log(first)-expected(1)) < 1e-7);
    assert(abs(log(last)-expected(end)) < 1e-7);

    % Dorduncu cikti icin de gizli son -3 dB kesisimi bulunmali.
    target = 10^(-3/20);
    expected3 = 10.4 + width*sqrt(log(2/(1+log(target))));
    assert(abs(log(last3)-expected3) < 1e-7);

    % Baslangicta hic grid kesisimi olmayan cukur ve tepe.
    for direction = [1, -1]
        Hfun = @(wi) exp(direction * ...
            (1-2*exp(-((log(wi)-3.4)/width).^2)));
        Hgrid = sampleResponse(w, Hfun);
        g = reshape(Hgrid, 1, []) - 1;
        assert(~any(g(1:end-1).*g(2:end) < 0));
        [~, ~, all] = svCrossings(w, Hgrid, Hfun, 0);
        assertRoots(all{1}, [3.4-offset, 3.4+offset]);
    end

    % Bilinen tek kesisim.
    Hfun = @(wi) exp(6-log(wi));
    [first, last, all] = svCrossings(w, sampleResponse(w,Hfun), Hfun);
    assertRoots(all{1}, 6);
    assert(abs(log(first)-6) < 1e-7 && abs(log(last)-6) < 1e-7);

    % Kesisim yoksa NaN ve bos liste; tek noktali grid de desteklenir.
    Hfun = @(wi) 0.01;
    for grid = {w, w(1)}
        [first, last, all, last3] = ...
            svCrossings(grid{1}, sampleResponse(grid{1},Hfun), Hfun);
        assert(isnan(first) && isnan(last) && isnan(last3));
        assert(isempty(all{1}));
    end

    % MIMO: ikinci tekil deger sabit, ilk dal gizli kesisimler iceriyor.
    Hfun = @(wi) diag([exp(profile(log(wi))), 0.01]);
    [first, last, all] = svCrossings(w, sampleResponse(w,Hfun), Hfun);
    assertRoots(all{1}, expected);
    assert(isnan(first(2)) && isnan(last(2)) && isempty(all{2}));

    % Parametrik hedef seviye ve -3 dB sonucunun yeniden kullanimi.
    Hfun = @(wi) target*exp(profile(log(wi)));
    [~, last, all, last3] = ...
        svCrossings(w, sampleResponse(w,Hfun), Hfun, -3);
    assertRoots(all{1}, expected);
    assert(last == last3);
    disp('Tum svCrossings kontrolleri gecti.');
end

function y = profile(x)
    y = max(-1, min(1, 6-x));
    for center = [1.4, 3.4]
        y = y - 2*exp(-((x-center)/0.32).^2);
    end
    for center = [8.4, 10.4]
        y = y + 2*exp(-((x-center)/0.32).^2);
    end
end

function Hgrid = sampleResponse(w, Hfun)
    H = Hfun(w(1));
    Hgrid = complex(zeros(size(H,1), size(H,2), numel(w)));
    for j = 1:numel(w)
        Hgrid(:,:,j) = Hfun(w(j));
    end
end

function assertRoots(actual, expected)
    assert(numel(actual) == numel(expected));
    assert(all(abs(log(actual(:))-expected(:)) < 1e-7));
end
