%% File Info.

%{

    simulate.m
    ----------
    % This code simulates the optimal life cycle decisions based on computed policies.

%}

%% Simulate class.

classdef simulate
    methods(Static)
        function sim = lc(par, sol)
            T = par.T;
            N = par.NN; % Number of agents
            a_sim = nan(T, N); % Asset simulation matrix
            c_sim = nan(T, N); % Consumption simulation matrix
            y_sim = nan(T, N); % Income simulation matrix
            alpha_sim = nan(T, N);
            r_sim = nan(T, N); % Stochastic interest rate
            
            % Initialize indices (ensure they are valid integers)
            a0_ind = round(rand(N,1) * (par.alen - 1)) + 1;
            t0_ind = ones(N,1); % Start from age 1
            y0_ind = round(rand(N,1) * (par.ylen - 1)) + 1;
            
            
            for i = 1:N
                for t = 1:T
                    c_sim(t, i) = sol.c(a0_ind(i), t0_ind(i), y0_ind(i)); % Get optimal consumption
                    a_sim(t, i) = sol.a(a0_ind(i), t0_ind(i), y0_ind(i)); % Get optimal assets
                    y_sim(t, i) = par.ygrid(y0_ind(i)); % Track income
                    r_t = par.r + 0.005 * randn();
                    r_sim(t, i) = r_t; % Track interest rate
                    alpha_sim(t, i) = sol.alpha(a0_ind(i), t0_ind(i), y0_ind(i));
                    
                    % Update indices for the next period
                    a0_ind(i) = round(a_sim(t, i) / max(par.agrid) * (par.alen - 1)) + 1;
                    y0_ind(i) = round(rand() * (par.ylen - 1)) + 1;
                end
            end
            
            sim.a = a_sim;
            sim.c = c_sim;
            sim.y = y_sim;
            sim.alpha = alpha_sim;
            sim.r = r_sim;
        end
    end
end