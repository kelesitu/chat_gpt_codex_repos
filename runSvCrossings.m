% Continuous-time MIMO example. Requires Control System Toolbox.
% Replace L with your own MIMO transfer function.
s = tf('s');
L = [10/(s+1), 0; 0, 2/(s+1)];

w = logspace(-3, 4, 1000);       % rad/s
Hgrid = freqresp(L, w);          % Reuse your existing response if available
Hfun = @(wi) evalfr(L, 1i*wi);   % Continuous-time models only
levelDB = -3;

[wFirst, wLast, wAll] = svCrossings(w, Hgrid, Hfun, levelDB);

disp(table((1:numel(wFirst))', wFirst, wLast, ...
    'VariableNames', {'SV', 'First_rad_s', 'Last_rad_s'}));

% Check this example against its analytical crossing frequencies.
expected = sqrt(([10; 2] / 10^(levelDB/20)).^2 - 1);
assert(all(abs(wFirst - expected)./expected < 1e-7));
assert(all(abs(wLast  - expected)./expected < 1e-7));

save('svSonuclari.mat', 'levelDB', 'wFirst', 'wLast', 'wAll');
