# Part 1: Linear Congruential Generator (LCG)
```matlab
% Generates the first 60 uniforms from the given LCG and prints the
% first three plus u51, u52, u53 (rounded only for display).
clear; clc;

x0 = 1000; a = 24693; c = 3517; K = 2^17; N = 60;
u = zeros(N,1); x = x0;
for i = 1:N
    x = mod(a*x + c, K);
    u(i) = x / K;
end

u_disp = round(u, 4); % round for printing only (do NOT use rounded values in sims)
fprintf('PART 1 — LCG CHECK\n');
fprintf('first3: %.4f, %.4f, %.4f\n', u_disp(1), u_disp(2), u_disp(3));
fprintf('u51=%.4f\nu52=%.4f\nu53=%.4f\n\n', u_disp(51), u_disp(52), u_disp(53));
```

# Part 2: Simulation
```matlab
% ---- LCG params (redeclare for clarity) ----
a = 24693; 
c = 3517; 
K = 2^17; 
LCG_x = 1000;              % seed

% ---- simulate 1000 customers ----
n = 1000;
W = zeros(n,1);

for i = 1:n
    W_i = 0;
    for attempt = 1:4
        % Dial
        W_i = W_i + 6;
        % One uniform for discrete outcome
        LCG_x = mod(a*LCG_x + c, K); 
        u2 = LCG_x / K;

        if u2 < 0.2
            % Busy -> voicemail 3s
            W_i = W_i + 3 + 1;      % + hang-up
        elseif u2 < 0.5
            % Unavailable -> rings 25s
            W_i = W_i + 25 + 1;     % + hang-up
        else
            % Available -> draw exponential with mean 12
            LCG_x = mod(a*LCG_x + c, K);
            u3 = LCG_x / K;
            x = -12 * log(1 - u3);  % inverse-CDF; guard not needed with this LCG

            if x <= 25
                W_i = W_i + x + 1;  % + hang-up, stop
                break
            else
                W_i = W_i + 25 + 1; % timed out, continue
            end
        end
    end
    W(i) = W_i;
end

% ---- Step 5 estimates ----
mu = mean(W);
if exist('quantile','file')
    q = quantile(W, [0.25 0.5 0.75]);
else
    q = prctile(W, [25 50 75]);   % fallback
end

p_le_15 = mean(W <= 15);
p_le_20 = mean(W <= 20);
p_le_30 = mean(W <= 30);
p_gt_40 = mean(W >  40);

w1 = 60; w2 = 80; w3 = 100;
p_gt_w1 = mean(W > w1);
p_gt_w2 = mean(W > w2);
p_gt_w3 = mean(W > w3);

t = [10 15 20 30 40 50 60 80 100 120 128];
F = arrayfun(@(tt) mean(W <= tt), t);

% ---- print results ----
fprintf('PART 2 — STEP 5 RESULTS\n');
fprintf('Mean: %.2f s\n', mu);
fprintf('Quartiles: Q1=%.2f s, Median=%.2f s, Q3=%.2f s\n', q(1), q(2), q(3));

fprintf('\nRequested probabilities:\n');
fprintf('P(W<=15) = %.3f\n', p_le_15);
fprintf('P(W<=20) = %.3f\n', p_le_20);
fprintf('P(W<=30) = %.3f\n', p_le_30);
fprintf('P(W>40)  = %.3f\n', p_gt_40);

fprintf('\nRight-tail points (w1=60, w2=80, w3=100):\n');
fprintf('P(W>60)  = %.3f\n', p_gt_w1);
fprintf('P(W>80)  = %.3f\n', p_gt_w2);
fprintf('P(W>100) = %.3f\n', p_gt_w3);

fprintf('\nEmpirical CDF checkpoints:\n');
for k = 1:numel(t)
    fprintf('t=%3d  F(W<=t)=%.3f\n', t(k), F(k));
end

```
