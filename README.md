# Part 1: Linear Congruential Generator (LCG)
```matlab
x0 = 1000; a = 24693; c = 3517; K = 2^17; N = 60;
u = zeros(N,1); x = x0;
for i = 1:N
 x = mod(a*x + c, K);
 u(i) = x / K;
end
% show first three u's (rounded)
fprintf('%.4f, %.4f, %.4f\n', u(1), u(2), u(3));
% show u51–u53 (rounded only)
fprintf('u51 = %.4f\n', u(51));
fprintf('u52 = %.4f\n', u(52));
fprintf('u53 = %.4f\n', u(53));
% discrete outcomes (print rounded u only)
for k = 51:53
    if u(k) < 0.2
        s = 'busy';
    elseif u(k) < 0.5
        s = 'unavailable';
    else
        s = 'available';
    end
    fprintf('Discrete: u%02d = %.4f -> %s\n', k, u(k), s);
end
% continuous times (print rounded u and rounded X only)
m = 12; epsSafe = 1e-12;
for k = 51:53
    uk = min(max(u(k), epsSafe), 1 - epsSafe);
    xk = -m * log(1 - uk); % 1 - uk (instead of log(uk)) bc some abt being more exact when uk near 0
    fprintf('Continuous: u%02d = %.4f -> X = %.4f s\n', k, u(k), xk);
end
```
