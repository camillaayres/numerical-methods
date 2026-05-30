function [price, se, CI_low, CI_high] = mc_uoc_antithetic(S0, K, r, sigma, T, U, n, M)
% mc_uoc_antithetic  Monte Carlo with antithetic variates for a discretely
%   monitored up-and-out barrier call option.
%
%   [price, se, CI_low, CI_high] = mc_uoc_antithetic(S0, K, r, sigma, T, U, n, M)
%
%   Inputs:
%       S0    - initial stock price
%       K     - strike price
%       r     - risk-free rate
%       sigma - volatility
%       T     - time to maturity
%       U     - upper barrier level
%       n     - number of monitoring dates
%       M     - total number of simulations (uses M/2 antithetic pairs)
%
%   Outputs:
%       price   - estimated option price E(C_n)
%       se      - standard error of the estimate
%       CI_low  - lower bound of 95% confidence interval
%       CI_high - upper bound of 95% confidence interval

% Time step and number of antithetic pairs
dt = T / n;
Mpairs = M / 2;

% Precompute constants for the GBM increments
drift = (r - 0.5*sigma^2) * dt;
vol   = sigma * sqrt(dt);
disc  = exp(-r*T);

% Preallocate pair-average payoffs
Y = zeros(Mpairs, 1);

for j = 1:Mpairs
    
    % One shock vector shared by both paths
    Z = randn(1, n);
    
    % Initialise both paths at S0
    S1 = S0;   S2 = S0;
    maxS1 = S0;  maxS2 = S0;
    ko1 = false; ko2 = false;
    
    for i = 1:n
        % Path 1: original shocks
        if ~ko1
            S1 = S1 * exp(drift + vol*Z(i));
            maxS1 = max(maxS1, S1);
            if maxS1 >= U, ko1 = true; end
        end
        % Path 2: antithetic shocks
        if ~ko2
            S2 = S2 * exp(drift - vol*Z(i));
            maxS2 = max(maxS2, S2);
            if maxS2 >= U, ko2 = true; end
        end
    end
    
    % Payoffs: zero if knocked out, discounted call otherwise
    if ko1, payoff1 = 0; else, payoff1 = disc * max(S1 - K, 0); end
    if ko2, payoff2 = 0; else, payoff2 = disc * max(S2 - K, 0); end
    
    % Antithetic pair average
    Y(j) = (payoff1 + payoff2) / 2;
    
end

% Price estimate, standard error, and 95% CI
price   = mean(Y);
se      = std(Y) / sqrt(Mpairs);
CI_low  = price - 1.96 * se;
CI_high = price + 1.96 * se;

end
