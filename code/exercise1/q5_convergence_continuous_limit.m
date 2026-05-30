% question5.m
% Computes E(LBn) for increasing n = 2^2, 2^3, ..., 2^10 and studies
% convergence to the continuous-monitoring limit E(LB_inf) from Thompson (1999).


clear; 
clc;

%Parameters 
S0    = 100;
r     = 0.04;
sigma = 0.3;
K     = 100;
T     = 1;
df    = exp(-r*T);

%Continuous limit E(LB_inf) 
% Formula (7): E(LB_inf) = S0/T * int_0^T e^{-r(T-t)} N(arg(t)) dt - K*e^{-rT}*N(-gamma/sqrt(T/3))
% where arg(t) = (-gamma + sigma*t - sigma*t^2/(2T)) / sqrt(T/3)
%       gamma  = (ln(K/S0) - (r-sigma^2/2)*T/2) / sigma

gamma = (log(K/S0) - (r - sigma^2/2)*T/2) / sigma;

integrand = @(t) exp(-r*(T - t)) .* normcdf((-gamma + sigma*t - sigma*t.^2/(2*T)) / sqrt(T/3));

LB_inf = (S0/T) * integral(integrand, 0, T) - K*df*normcdf(-gamma/sqrt(T/3));

fprintf('E(LB_inf) = %.8f\n\n', LB_inf);

%Discrete E(LBn) for n = 2^2 to 2^10
powers = 2:10;
n_vec  = 2.^powers;
LB_n   = zeros(size(n_vec));

for idx = 1:length(n_vec)
    n  = n_vec(idx);
    dt = T / n;
    k  = (1:n)';

    % Shared quantities from equation (4)
    sigma_bar = sigma * sqrt((2*n + 1) / (3*n));
    T_bar     = (n + 1) * dt / 2;
    v         = sigma_bar * sqrt(T_bar);

    % Per-date quantities
    mu_k     = (r - 0.5*sigma^2) * k * dt;
    sig_k    = sigma * sqrt(k * dt);
    % Correct a_k: linear in k*(n+1-(k+1)/2), not sqrt
    a_k      = sigma * sqrt(dt) .* k .* (n + 1 - (k+1)/2) / sqrt(n*(n+1)*(2*n+1)/6);
    exp_term = exp(mu_k + 0.5*sig_k.^2);

    b = (log(S0/K) + (r - 0.5*sigma^2)*T_bar) / v;

    % Curran formula (3)
    LB_n(idx) = (S0*df/n) * sum(exp_term .* normcdf(b + a_k)) - K*df*normcdf(b);
end

%Print convergence table 
fprintf('%-6s | %-14s | %-18s | %-8s\n', 'n', 'E(LBn)', '|E(LBn)-E(LBinf)|', 'Ratio');
fprintf('%s\n', repmat('-', 1, 55));
prev_err = NaN;
for idx = 1:length(n_vec)
    err = abs(LB_n(idx) - LB_inf);
    if isnan(prev_err)
        ratio_str = '---';
    else
        ratio_str = sprintf('%.3f', prev_err/err);
    end
    fprintf('%-6d | %-14.8f | %-18.8f | %s\n', n_vec(idx), LB_n(idx), err, ratio_str);
    prev_err = err;
end

% Plot
fig = figure;
set(fig, 'Units', 'centimeters', 'Position', [0 0 18 7]); 

subplot(1,2,1);
plot(powers, LB_n, 'b-o', 'LineWidth', 1.5, 'MarkerFaceColor', 'b', 'MarkerSize', 5);
hold on;
yline(LB_inf, 'r--', 'LineWidth', 1.5);
xlabel('log_2(n)'); ylabel('Value');
title('E(LB_n) convergence to E(LB_\infty)');
legend('E(LB_n)', 'E(LB_\infty)', 'Location', 'northeast');
set(gca, 'XTick', powers); grid on;

subplot(1,2,2);
semilogy(powers, abs(LB_n - LB_inf), 'b-o', 'LineWidth', 1.5, ...
         'MarkerFaceColor', 'b', 'MarkerSize', 5);
xlabel('log_2(n)'); ylabel('|E(LB_n) - E(LB_\infty)|');
title('Absolute error (log scale)');
set(gca, 'XTick', powers); grid on;
exportgraphics(fig, 'fig_q5_convergence.pdf', 'ContentType', 'vector');
