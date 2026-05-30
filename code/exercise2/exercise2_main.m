% exercise2_main.m
% Main script for Exercise 2: Discretely monitored up-and-out barrier call
% option pricing via Monte Carlo simulation with antithetic variates.
%
% This script:
%   1) Computes the continuous-monitoring closed-form benchmark
%   2) Runs antithetic MC for all (U, n) pairs
%   3) Runs crude MC for selected n to quantify the variance reduction
%   4) Produces result tables and convergence figure
%
% Coursework parameters:
%   S0 = 110, K = 100, r = 0.05, sigma = 0.1, T = 2
%   U = [160 170], n = 2^2, 2^3, ..., M = 10^6

clear;

%% Parameters
S0    = 110;
K     = 100;
r     = 0.05;
sigma = 0.1;
T     = 2;
M     = 1e6;
Uvec  = [160 170];

% Monitoring dates: n = 2^2, ..., 2^10
powers = 2:10;
nvec   = 2.^powers;

%% Continuous-monitoring benchmarks
fprintf('=== Continuous-monitoring benchmark (formula (8)) ===\n');
f_cont = zeros(length(Uvec), 1);
for iu = 1:length(Uvec)
    f_cont(iu) = fUOC_continuous(T, S0, K, Uvec(iu), sigma, r);
    fprintf('  U = %d:  f_UOC = %.6f\n', Uvec(iu), f_cont(iu));
end

%% Antithetic MC for all (U, n) pairs
nU = length(Uvec);
nN = length(nvec);
prices  = zeros(nU, nN);
SEs     = zeros(nU, nN);
CI_lows = zeros(nU, nN);
CI_his  = zeros(nU, nN);

for iu = 1:nU
    U = Uvec(iu);
    fprintf('\n=== Antithetic MC: U = %d ===\n', U);
    for jn = 1:nN
        n = nvec(jn);
        tic;
        [p, se, cl, ch] = mc_uoc_antithetic(S0, K, r, sigma, T, U, n, M);
        elapsed = toc;
        
        prices(iu, jn)  = p;
        SEs(iu, jn)     = se;
        CI_lows(iu, jn) = cl;
        CI_his(iu, jn)  = ch;
        
        fprintf('  n=%4d | Price=%.6f | SE=%.6f | CI=[%.6f,%.6f] | %.1fs\n', ...
                n, p, se, cl, ch, elapsed);
    end
end

%% Crude MC comparison for selected n values
compare_n = [4, 64];
fprintf('\n=== Crude MC comparison ===\n');
fprintf('%4s %4s %12s %10s %12s %10s %8s\n', ...
        'U', 'n', 'Crude SE', 'AV SE', 'Crude CIw', 'AV CIw', 'VarRatio');

for iu = 1:nU
    U = Uvec(iu);
    for cn = compare_n
        jn = find(nvec == cn);
        
        tic;
        [~, se_cr, cl_cr, ch_cr] = mc_uoc_crude(S0, K, r, sigma, T, U, cn, M);
        toc;
        
        se_av  = SEs(iu, jn);
        ciw_cr = ch_cr - cl_cr;
        ciw_av = CI_his(iu, jn) - CI_lows(iu, jn);
        vr     = (se_cr / se_av)^2;
        
        fprintf('%4d %4d %12.6f %10.6f %12.6f %10.6f %8.2f\n', ...
                U, cn, se_cr, se_av, ciw_cr, ciw_av, vr);
    end
end

%% Print result tables
for iu = 1:nU
    U = Uvec(iu);
    fprintf('\n--- Table: U = %d (benchmark = %.6f) ---\n', U, f_cont(iu));
    fprintf('%6s %12s %10s %12s %12s %10s %10s\n', ...
            'n', 'MC Price', 'Std Error', 'CI Lower', 'CI Upper', 'CI Width', 'Gap');
    for jn = 1:nN
        ci_w = CI_his(iu,jn) - CI_lows(iu,jn);
        gap  = prices(iu,jn) - f_cont(iu);
        fprintf('%6d %12.6f %10.6f %12.6f %12.6f %10.6f %10.6f\n', ...
                nvec(jn), prices(iu,jn), SEs(iu,jn), ...
                CI_lows(iu,jn), CI_his(iu,jn), ci_w, gap);
    end
end

%% Convergence figure
figure;
hold on;
markers = {'o-', 's-'};
colours = {'b', 'r'};
for iu = 1:nU
    errorbar(powers, prices(iu,:), ...
             prices(iu,:) - CI_lows(iu,:), CI_his(iu,:) - prices(iu,:), ...
             markers{iu}, 'Color', colours{iu}, 'LineWidth', 1.4, ...
             'MarkerFaceColor', colours{iu}, 'MarkerSize', 5, 'CapSize', 3);
end
for iu = 1:nU
    yline(f_cont(iu), '--', sprintf('f_{UOC} (U=%d)', Uvec(iu)), ...
           'Color', colours{iu}, 'LineWidth', 1, ...
           'LabelHorizontalAlignment', 'left');
end
xlabel('log_2(n)');
ylabel('Estimated option price');
title('Convergence of discrete UOC price to continuous benchmark');
legend(sprintf('U = %d (discrete MC)', Uvec(1)), ...
       sprintf('U = %d (discrete MC)', Uvec(2)), 'Location', 'northeast');
set(gca, 'XTick', powers);
grid on; hold off;
saveas(gcf, 'exercise2_convergence.png');
fprintf('\nFigure saved. Exercise 2 complete.\n');
