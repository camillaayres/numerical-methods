% File: question4.m
% Purpose:
%   Monte Carlo estimation of arithmetic Asian call option prices AND LR deltas
%   using two control variates:
%     (a) Control variate = LB_n (Curran lower bound)
%     (b) Control variate = G_n  (geometric Asian, Kemna-Vorst)
%
%   For each (K,n) pair, reports: estimate, standard error, 95% CI, and
%   efficiency ratio E(K,n) = (t_a * SE_a^2) / (t_b * SE_b^2).
%   E < 1 means method (a) is more efficient than method (b).
%
%   Closed-form expectations from Questions 2 and 3 are used as control means.

clear; clc;

%% --- Parameters ---
S0    = 100;
r     = 0.04;
sigma = 0.3;
T     = 1;

K_vec = [90 100 110];
n_vec = [4 12 50];

M  = 1e5;       % number of Monte Carlo paths
z  = 1.96;      % 95% CI z-score
df = exp(-r*T); % discount factor

rng(12345, 'twister'); % fix seed for reproducibility

%% --- Storage: one row per (K,n) pair ---
% Columns: [K, n, est_a, se_a, ciL_a, ciU_a, t_a, est_b, se_b, ciL_b, ciU_b, t_b, Eff]
nRows      = numel(K_vec) * numel(n_vec);
price_rows = zeros(nRows, 13);
delta_rows = zeros(nRows, 13);
ix = 0;

%% --- Main loop ---
for ni = 1:numel(n_vec)
    n  = n_vec(ni);
    dt = T / n;
    k  = (1:n)'; % column vector of monitoring indices

    % Shared quantities from equation (4)
    sigma_bar = sigma * sqrt((2*n + 1) / (3*n));
    T_bar     = (n + 1) * dt / 2;
    v         = sigma_bar * sqrt(T_bar);

    % Per-date Curran quantities
    mu_k     = (r - 0.5*sigma^2) .* (k * dt);
    sig_k    = sigma .* sqrt(k * dt);
    % a_k from Curran (1994): Cov(log S_k, log G_n) / sqrt(Var(log G_n))
    % = sigma*sqrt(dt) * k*(n+1-(k+1)/2) / sqrt(n*(n+1)*(2n+1)/6)
    % Note: numerator is LINEAR in k*(n+1-(k+1)/2), not sqrt
    a_k      = sigma * sqrt(dt) .* (k .* (n + 1 - (k+1)/2)) ./ sqrt(n*(n+1)*(2*n+1)/6);
    exp_term = exp(mu_k + 0.5 * sig_k.^2); % used in Curran formulas

    for ki = 1:numel(K_vec)
        K  = K_vec(ki);
        ix = ix + 1;

        % --- Simulate M GBM paths ---
        Z    = randn(M, n);
        logS = log(S0) + cumsum((r - 0.5*sigma^2)*dt + sigma*sqrt(dt).*Z, 2);
        S    = exp(logS);

        Amean = mean(S, 2);         % arithmetic average per path
        Gmean = exp(mean(logS, 2)); % geometric average per path

        % --- Discounted payoff samples ---
        A  = df .* max(Amean - K, 0);           % arithmetic Asian call, eq. (1)
        LB = df .* (Amean - K) .* (Gmean > K);  % lower bound payoff,    eq. (2)
        G  = df .* max(Gmean - K, 0);           % geometric Asian call,   eq. (5)

        % --- Closed-form expectations from Questions 2-3 ---
        b = (log(S0/K) + (r - 0.5*sigma^2)*T_bar) / v;              % eq. (4)
        d = (log(S0/K) + (r - 0.5*sigma^2 + sigma_bar^2)*T_bar) / v; % eq. (4)

        % E[LBn]: Curran formula (3)
        ELB = (S0*df/n) * sum(exp_term .* normcdf(b + a_k)) - K*df*normcdf(b);

        % E[Gn]: Kemna-Vorst formula (6)
        alpha_KV = exp((r - 0.5*sigma^2 + 0.5*sigma_bar^2)*T_bar - r*T);
        EG       = S0*alpha_KV*normcdf(d) - K*df*normcdf(d - v);

        % dE[LBn]/dS0: Question 2
        part1   = (df/n) * sum(exp_term .* normcdf(b + a_k));
        part2   = (df/(S0*v)) * ((S0/n)*sum(exp_term .* normpdf(b + a_k)) - K*normpdf(b));
        ELBSENS = part1 + part2;

        % dE[Gn]/dS0: Question 2
        EGSENS = alpha_KV*normcdf(d) + (1/v)*(alpha_KV*normpdf(d) - (K/S0)*df*normpdf(d - v));

        % --- LR delta score ---
        S1    = S(:, 1);
        zeta  = (log(S1/S0) - (r - 0.5*sigma^2)*dt) / (sigma*sqrt(dt));
        score = zeta / (S0 * sigma * sqrt(dt));

        % LR delta samples = payoff * score
        dA  = A  .* score;
        dLB = LB .* score;
        dG  = G  .* score;

        % --- CV estimation: prices ---
        tStart = tic;
        [pLB, seLB, ciLB] = cv_estimator(A, LB, ELB, z);
        tLB_p = toc(tStart);

        tStart = tic;
        [pG, seG, ciG] = cv_estimator(A, G, EG, z);
        tG_p = toc(tStart);

        % --- CV estimation: deltas ---
        tStart = tic;
        [dLB_est, se_dLB, ci_dLB] = cv_estimator(dA, dLB, ELBSENS, z);
        tLB_d = toc(tStart);

        tStart = tic;
        [dG_est, se_dG, ci_dG] = cv_estimator(dA, dG, EGSENS, z);
        tG_d = toc(tStart);

        % --- Efficiency ratios ---
        Eff_p = (tLB_p * seLB^2)   / (tG_p * seG^2);
        Eff_d = (tLB_d * se_dLB^2) / (tG_d * se_dG^2);

        % --- Store ---
        price_rows(ix,:) = [K, n, pLB, seLB, ciLB(1), ciLB(2), tLB_p, ...
                                   pG,  seG,  ciG(1),  ciG(2),  tG_p,  Eff_p];
        delta_rows(ix,:) = [K, n, dLB_est, se_dLB, ci_dLB(1), ci_dLB(2), tLB_d, ...
                                   dG_est,  se_dG,  ci_dG(1),  ci_dG(2),  tG_d,  Eff_d];
    end
end

%% --- Display results ---
cols = {'K','n','Estimate_CV_LB','SE','CI_Lower','CI_Upper','EfficiencyRatio'};

PriceTable_LB = array2table(price_rows(:,[1 2 3 4 5 6 13]),  'VariableNames', cols);
PriceTable_G  = array2table(price_rows(:,[1 2 8 9 10 11 13]), 'VariableNames', cols);
DeltaTable_LB = array2table(delta_rows(:,[1 2 3 4 5 6 13]),  'VariableNames', cols);
DeltaTable_G  = array2table(delta_rows(:,[1 2 8 9 10 11 13]), 'VariableNames', cols);

disp('=== Prices: control variate (a) = LB_n ===');      disp(PriceTable_LB);
disp('=== Prices: control variate (b) = G_n ===');       disp(PriceTable_G);
disp('=== LR Deltas: control variate (a) = LB_n delta ==='); disp(DeltaTable_LB);
disp('=== LR Deltas: control variate (b) = G_n delta ===');  disp(DeltaTable_G);

%% ===== Local function =====
function [mu, se, ci] = cv_estimator(Y, X, EX, z)
% CV estimator: Y_cv = Y - beta*(X - E[X]),  beta = Cov(Y,X)/Var(X)
    Y    = Y(:);  X = X(:);
    M    = numel(Y);
    Xc   = X - mean(X);
    Yc   = Y - mean(Y);
    beta = (Xc'*Yc) / (Xc'*Xc + eps); % eps guards against zero variance
    Ycv  = Y - beta*(X - EX);
    mu   = mean(Ycv);
    se   = std(Ycv, 0) / sqrt(M);     % std with ddof=1
    ci   = [mu - z*se, mu + z*se];
end
