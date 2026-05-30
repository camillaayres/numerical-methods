function price = fUOC_continuous(T, S0, K, U, sigma, r)
% fUOC_continuous  Closed-form price of a continuously monitored
%   up-and-out barrier call option under GBM (Shreve, 2004).
%
%   price = fUOC_continuous(T, S0, K, U, sigma, r)
%
%   Inputs:
%       T     - time to maturity
%       S0    - initial stock price
%       K     - strike price
%       U     - upper barrier level
%       sigma - volatility
%       r     - risk-free interest rate
%
%   Output:
%       price - closed-form UOC option price (formula (8) in coursework)

% d_+(x) and d_-(x) as defined in equation (9)
dp = @(x) (log(x) + (r + 0.5*sigma^2)*T) / (sigma*sqrt(T));
dm = @(x) (log(x) + (r - 0.5*sigma^2)*T) / (sigma*sqrt(T));

% Exponent appearing in the reflection terms
alpha = -2*r / sigma^2;

% Assemble the four terms of formula (8)
term1 = S0 * ( normcdf(dp(S0/K)) - normcdf(dp(S0/U)) );
term2 = -K * exp(-r*T) * ( normcdf(dm(S0/K)) - normcdf(dm(S0/U)) );
term3 = -U * (S0/U)^alpha ...
        * ( normcdf(dp(U^2/(K*S0))) - normcdf(dp(U/S0)) );
term4 =  K * exp(-r*T) * (S0/U)^(1 + alpha) ...
        * ( normcdf(dm(U^2/(K*S0))) - normcdf(dm(U/S0)) );

price = term1 + term2 + term3 + term4;

end
