function [price, se, CI_low, CI_high] = mc_uoc_crude(S0, K, r, sigma, T, U, n, M)
% mc_uoc_crude  Crude (standard) Monte Carlo for a discretely monitored
%   up-and-out barrier call option. No variance reduction is applied.
%
%   [price, se, CI_low, CI_high] = mc_uoc_crude(S0, K, r, sigma, T, U, n, M)
%
%   Used for comparison against the antithetic variate estimator.

dt = T / n;
drift = (r - 0.5*sigma^2) * dt;
vol   = sigma * sqrt(dt);
disc  = exp(-r*T);

% Preallocate individual payoffs
C = zeros(M, 1);

for j = 1:M
    
    Z = randn(1, n);
    S = S0;
    maxS = S0;
    ko = false;
    
    for i = 1:n
        if ~ko
            S = S * exp(drift + vol*Z(i));
            maxS = max(maxS, S);
            if maxS >= U, ko = true; end
        end
    end
    
    if ko
        C(j) = 0;
    else
        C(j) = disc * max(S - K, 0);
    end
    
end

price   = mean(C);
se      = std(C) / sqrt(M);
CI_low  = price - 1.96 * se;
CI_high = price + 1.96 * se;

end
