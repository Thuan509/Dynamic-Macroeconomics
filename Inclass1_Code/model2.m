%% File Info.

%{
    model.m
    -------
    This code sets up the Life Cycle Model parameters with stochastic income during working years.
%}

%% Model class.

classdef model2
    methods(Static)
        %% Set up structure array for model parameters and set the simulation parameters.
        
        function par = setup()            
            %% Structure array for model parameters.
            
            par = struct();
            
            %% Preferences.
            
            par.beta = 0.96; % Discount factor
            par.sigma = 2.00; % CRRA risk aversion parameter
            
            assert(par.beta > 0 && par.beta < 1.00, 'Discount factor should be between 0 and 1.\n')
            assert(par.sigma > 0, 'CRRA should be at least 0.\n')

            %% Income and Technology Parameters

            par.ybar = 1.0; % Base income level (for stochastic process)
            par.kappa = 0.3; % Pension fraction of income in last working period
            par.r = 0.1; % Fixed interest rate
            par.a0 = 5; % Initial wealth endowment (a0 > 0)
            par.aT = 0; % Terminal condition (no bequests)

            % Stochastic income process parameters
            par.rho = 0.9; % Persistence parameter for AR(1) process (rho < 1)
            par.sigma_e = 0.2; % Standard deviation of income shocks
            
            assert(par.ybar > 0, 'The base income should be larger than 0.\n')
            assert(par.kappa > 0 && par.kappa < 1, 'The pension fraction should be between 0 and 1.\n')
            assert(par.a0 > 0, 'The initial endowment should be greater than 0.\n')
            assert(par.aT == 0, 'The terminal wealth should be equal to 0.\n')
            assert(par.rho < 1, 'The persistence parameter should be less than 1.\n')
            assert(par.sigma_e > 0, 'The standard deviation of shocks should be positive.\n')

            %% Simulation parameters.

            % Option 1: Shorter life (for testing)
            par.T = 4; % Total lifespan (4 periods)
            par.tr = 2; % Retirement begins at period 2 (work periods: 0,1; retirement periods: 2,3)
            
            % Option 2: Longer life (uncomment to use)
            % par.T = 60; % Total lifespan (60 periods)
            % par.tr = 41; % Retirement begins at period 41 (work periods: 0-40; retirement: 41-59)
            
            par.seed = 2025; % Random seed for simulation
        end
        
        %% Generate state grids.
        
        function par = gen_grids(par)
            %% Wealth grid.
            par.alen = 300; % Grid size for a.
            par.amax = 30; % Upper bound for a.
            par.amin = 0; % Minimum a.
            
            assert(par.alen > 5, 'Grid size for a should be positive and greater than 5.\n')
            assert(par.amax > par.amin, 'Minimum a should be less than maximum value.\n')
            
            par.agrid = linspace(par.amin, par.amax, par.alen)'; % Equally spaced grid for a.
            
            %% Income grid (discretize AR(1) using Tauchen method)
            par.ylen = 7; % Number of income states
            par.m = 3; % Number of standard deviations for income grid
            
            % Use Tauchen method to discretize income process
            [par.ygrid, par.ytrans] = model2.tauchen(par.ylen, par.rho, par.sigma_e, par.m);
            par.ygrid = exp(par.ygrid); % Convert from log to level
            
            % Scale income to have mean ybar
            mean_y = sum(par.ygrid .* sum(par.ytrans, 2) / par.ylen);
            par.ygrid = par.ygrid * (par.ybar / mean_y);
            
            fprintf('Income grid created with %d states\n', par.ylen)
            fprintf('Income states: [%.2f, %.2f, ..., %.2f]\n', par.ygrid(1), par.ygrid(2), par.ygrid(end))
        end
        
        %% Tauchen method for discretizing AR(1) process
        function [zgrid, P] = tauchen(n, rho, sigma, m)
            % Tauchen method for discretizing AR(1) process
            % log(y_t) = rho*log(y_{t-1}) + e_t, e_t ~ N(0, sigma^2)
            
            % Parameters
            sigma_z = sigma / sqrt(1 - rho^2); % Unconditional standard deviation
            
            % Create grid for z
            z_max = m * sigma_z;
            z_min = -z_max;
            zgrid = linspace(z_min, z_max, n)';
            w = zgrid(2) - zgrid(1); % Width of each interval
            
            % Compute transition matrix
            P = zeros(n, n);
            
            % For each current state
            for i = 1:n
                % For each possible next state
                for j = 1:n
                    if j == 1 % Lowest state
                        P(i, j) = normcdf((zgrid(j) + w/2 - rho*zgrid(i)) / sigma);
                    elseif j == n % Highest state
                        P(i, j) = 1 - normcdf((zgrid(j) - w/2 - rho*zgrid(i)) / sigma);
                    else % Middle states
                        P(i, j) = normcdf((zgrid(j) + w/2 - rho*zgrid(i)) / sigma) - ...
                                  normcdf((zgrid(j) - w/2 - rho*zgrid(i)) / sigma);
                    end
                end
            end
            
            % Ensure rows sum to 1 (due to numerical issues)
            for i = 1:n
                P(i,:) = P(i,:) / sum(P(i,:));
            end
        end
        
        %% Utility function.
        
        function u = utility(c, par)
            %% CRRA utility.
            
            if par.sigma == 1
                u = log(c); % Log utility.
            else
                u = (c.^(1-par.sigma))./(1-par.sigma); % CRRA utility.
            end
        end
        %% Add these methods to the model class in model.m

% Add these methods to the model class:

        %% Interpolation functions for simulation
        
        function [idx_low, weight] = interp_index(x, grid)
            % Find interpolation index and weight for value x in grid
            
            n = length(grid);
            
            if x <= grid(1)
                idx_low = 1;
                weight = 1;
            elseif x >= grid(n)
                idx_low = n-1;
                weight = 0;
            else
                idx_low = find(grid <= x, 1, 'last');
                weight = (grid(idx_low+1) - x) / (grid(idx_low+1) - grid(idx_low));
            end
        end
        
        function [c, a_next] = interp_policy(a, y_idx, t, agrid, a_pol, c_pol)
            % Interpolate policy functions for a given state
            
            % Find position in asset grid
            [a_idx_low, weight] = model2.interp_index(a, agrid);
            a_idx_high = min(a_idx_low + 1, length(agrid));
            
            % Interpolate consumption
            c_low = c_pol(a_idx_low, y_idx, t);
            c_high = c_pol(a_idx_high, y_idx, t);
            c = weight * c_low + (1 - weight) * c_high;
            
            % Interpolate next period assets
            a_next_low = a_pol(a_idx_low, y_idx, t);
            a_next_high = a_pol(a_idx_high, y_idx, t);
            a_next = weight * a_next_low + (1 - weight) * a_next_high;
        end
        
    end
end