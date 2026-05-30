% SMM313 Group Coursework - Exercise 1, Task 3


% Computes exact closed-form values of:
% a) Expected lower bound E(LBn) and its sensitivity w.r.t. S0 (Curran (3))
% b) Geometric Asian call price E(Gn) and its sensitivity w.r.t. S0 (Kemna-Vorst (6)) for all (K, n) pairs.

clear
clc
clear all


% Parameters
S0    = 100;
r     = 0.04;
sigma = 0.3;
T     = 1;
K_vec = [90, 100, 110];
n_vec = [4, 12, 50];
 

% Loop over all (K, n) pairs 
fprintf('%-6s %-6s | %-12s %-12s | %-12s %-12s\n','K', 'n', 'E(LBn)', 'dE(LBn)/dS0', 'E(Gn)', 'dE(Gn)/dS0');
fprintf('%s\n', repmat('-', 1, 70));

for ki = 1:length(K_vec)
    K = K_vec(ki);
    for ni = 1:length(n_vec)
        n = n_vec(ni);
        dt = T / n; 

        % Shared quantities from equation (4)
        sigma_bar = sigma * sqrt((2*n + 1) / (3*n));
        T_bar     = (n + 1) * dt / 2;
        b         = (log(S0/K) + (r - sigma^2/2) * T_bar) / (sigma_bar * sqrt(T_bar));

        % Per monitoring date quantities
        mu_k    = zeros(1, n);
        sigma_k = zeros(1, n);
        a_k     = zeros(1, n);
        for k = 1:n
            mu_k(k)    = (r - sigma^2/2) * k * dt;
            sigma_k(k) = sigma * sqrt(k * dt);
           a_k(k) = sigma * sqrt(dt) * k * (n + 1 - (k+1)/2) / sqrt(n*(n+1)*(2*n+1)/6);
        end

        % a) Expected lower bound E(LBn) -- Curran formula (3)
        sum_LB      = sum(exp(mu_k + sigma_k.^2 / 2) .* normcdf(b + a_k));
        ELBn        = (S0 * exp(-r*T) / n) * sum_LB - K * exp(-r*T) * normcdf(b);




        % Sensitivity dE(LBn)/dS0
      

        %   Part 1: (e^{-rT}/n) * sum_k exp(mu_k + sig_k^2/2) * N(b + a_k)
        %           = ELBn/S0 + K*e^{-rT}*N(b)/S0   [equiv. the sum term / S0]
        

        %   Part 2: (e^{-rT} / (S0 * sigma_bar * sqrt(T_bar))) *
        %           [ (S0/n)*sum_k exp(mu_k+sig_k^2/2)*phi(b+a_k) - K*phi(b) ]


        db_dS0      = 1 / (S0 * sigma_bar * sqrt(T_bar));
        sum_phi_LB  = sum(exp(mu_k + sigma_k.^2 / 2) .* normpdf(b + a_k));

        part1       = (exp(-r*T) / n) * sum_LB;                         % d/dS0 of explicit S0 factor
        part2_inner = (S0 / n) * sum_phi_LB - K * normpdf(b);
        part2       = exp(-r*T) * db_dS0 * part2_inner;            % d/dS0 from b inside N(.)
        dELBn_dS0   = part1 + part2;


        % b) Geometric Asian option price E(Gn) -- Kemna-Vorst (6)
      
        d    = (log(S0/K) + (r - sigma^2/2 + sigma_bar^2) * T_bar) / (sigma_bar * sqrt(T_bar));
        alpha_KV = exp((r - sigma^2/2 + sigma_bar^2/2) * T_bar - r*T); % constant w.r.t. S0
        EGn  = S0 * alpha_KV * normcdf(d) - K * exp(-r*T) * normcdf(d - sigma_bar*sqrt(T_bar));

        % Sensitivity dE(Gn)/dS0 


        dd_dS0     = 1 / (S0 * sigma_bar * sqrt(T_bar));
        part1_G    = alpha_KV * normcdf(d);
        part2_G    = (1 / (sigma_bar * sqrt(T_bar))) * ...
                     (alpha_KV * normpdf(d) - (K/S0) * exp(-r*T) * normpdf(d - sigma_bar*sqrt(T_bar)));
        dEGn_dS0   = part1_G + part2_G;

        
        % Print results
        fprintf('%-6d %-6d | %-12.6f %-12.6f | %-12.6f %-12.6f\n',K, n, ELBn, dELBn_dS0, EGn, dEGn_dS0);

        % Storing results
        results(ki, ni).K          = K;
        results(ki, ni).n          = n;
        results(ki, ni).ELBn       = ELBn;
        results(ki, ni).dELBn_dS0  = dELBn_dS0;
        results(ki, ni).EGn        = EGn;
        results(ki, ni).dEGn_dS0   = dEGn_dS0;
        results(ki, ni).b          = b;
        results(ki, ni).d          = d;
        results(ki, ni).sigma_bar  = sigma_bar;
        results(ki, ni).T_bar      = T_bar;
        results(ki, ni).mu_k       = mu_k;
        results(ki, ni).sigma_k    = sigma_k;
        results(ki, ni).a_k        = a_k;
    end
end