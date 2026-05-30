% fig_pw_discontinuity.m
% Plots the LB_n payoff as a function of S0 for a fixed simulated path,
% illustrating the jump discontinuity at S0* where Gmean = K.

clear; clc;

K  = 100;  r = 0.04;  sigma = 0.3;  T = 1;  n = 4;
df = exp(-r*T);

% Fix log-increments so the boundary falls exactly at S0* = K = 100.
% Condition: mean(cumsum(log_inc)) = 0
log_inc = [0.3, -0.4, 0.2, -0.4];
cs      = cumsum(log_inc);

S0_vec    = linspace(88, 113, 600);
LB_payoff = zeros(size(S0_vec));

for i = 1:length(S0_vec)
    S0    = S0_vec(i);
    S     = S0 * exp(cs);
    Amean = mean(S);
    Gmean = exp(mean(log(S)));
    LB_payoff(i) = df * (Amean - K) * (Gmean > K);
end

S0_star = K;  % boundary by construction

fig = figure;
set(fig, 'Units', 'centimeters', 'Position', [0 0 10 6]);

plot(S0_vec, LB_payoff, 'r-', 'LineWidth', 2);
hold on;
xline(S0_star, ':', 'Color', [0.4 0.4 0.4], 'LineWidth', 1.4);
xlabel('$S_0$', 'Interpreter', 'latex', 'FontSize', 12);
ylabel('Discounted payoff', 'FontSize', 11);
title('$\mathrm{LB}_n$ payoff vs.\ $S_0$ (path fixed)', ...
      'Interpreter', 'latex', 'FontSize', 11);
legend({'$\mathrm{LB}_n$ payoff', 'Boundary $S_0^*$'}, ...
       'Interpreter', 'latex', 'FontSize', 10, 'Location', 'northwest');
grid on;
xlim([88 113]);

exportgraphics(fig, 'fig_pw_discontinuity.pdf', 'ContentType', 'vector');
