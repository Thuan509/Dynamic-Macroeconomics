%% File Info.
%{
    model2.m
    -------
    This code sets up the model.
%}

%% Model class.
classdef model2
    methods(Static)
        %% Set up structure array for model parameters and set the simulation parameters.
        function par = setup()            
            %% Structure array for model parameters.
            par = struct();
            
            %% Preferences.
            par.T = 61; % Last period of life.
            par.tr = 41; % First period of retirement.
            
            par.beta = 0.94; % Discount factor.
            par.sigma = 2.0; % CRRA risk aversion parameter.
            par.gamma = 1.00; % Weight on leisure: Higher values mean that leisure has a higher weight in the utility function.
            par.nu = 0.5; % Frisch Elasticity: Higher values of this mean that the labor choice becomes more sensitive to productivity shocks
            
            assert(par.T > par.tr, 'Cannot retire after dying.\n');
            assert(par.beta > 0.0 && par.beta < 1.0, 'Discount factor should be between 0 and 1.\n');
            assert(par.gamma > 0.0, 'CRRA should be positive.\n');
            assert(par.nu > 0,'The sensitivity to productivity shocks should be positive.\n')
            assert(par.gamma > 0,'The weight on leisure should be at least 0.\n')

            %% Prices and Income.
            par.r = 0.03; % Interest rate.
            par.ybar = 10.0; % Exogenous income.
            par.kappa = 0.6; % Share of income as pension.

            %% Load G_t values from CSV
            gt_data = readtable(fullfile(pwd, 'Gt_values.csv')); % Use pwd instead of main
            par.age_groups = gt_data.Age; % Store age indices
            par.Gt = gt_data.G_t; % Store corresponding G_t values
            par.Gmat = par.Gt;
            par.Gmat = par.Gmat / par.Gmat(1); % Normalize by first period's value
            par.Gmat = par.Gmat(1:61); % Only take values from 1 to 50 (i.e., )

            par.sigma_eps = 0.07; % Std. dev of productivity shocks.
            par.rho = 0.85; % Persistence of AR(1) process.
            par.mu = 0.0; % Intercept of AR(1) process.

            assert(par.ybar > 0.0, 'Income must be positive.\n');
            assert(par.kappa >= 0.0 && par.kappa <= 1.0, 'Pension share should be from 0 to 1.\n');
            assert(par.sigma_eps > 0, 'The standard deviation of the shock must be positive.\n');
            assert(abs(par.rho) < 1, 'The persistence must be less than 1 for stationarity.\n');

            %% Simulation parameters.
            par.seed = 2025; % Seed for simulation.
            par.TT = 61; % Number of time periods.
            par.NN = 10; % Number of people.
        end
        
        %% Generate state grids.
        function par = gen_grids(par)
            %% Capital grid.
            par.alen = 50; % Grid size for a.
            par.amax = 30.0; % Upper bound for a.
            par.amin = 0.0; % Minimum a.
        
            assert(par.alen > 5, 'Grid size for a should be positive and greater than 5.\n');
            assert(par.amax > par.amin, 'Minimum a should be less than maximum value.\n');
        
            par.agrid = linspace(par.amin, par.amax, par.alen)'; % Equally spaced, linear grid.
        
            %% Discretized income process.
            par.ylen = 7; % Grid size for y.
            par.m = 3; % Scaling parameter for Tauchen.
        
            assert(par.ylen > 3, 'Grid size for y should be positive and greater than 3.\n');
            assert(par.m > 0, 'Scaling parameter for Tauchen should be positive.\n');
        
            [ygrid, pmat] = model.tauchen(par.mu, par.rho, par.sigma_eps, par.ylen, par.m); % Tauchen's Method
            par.ygrid = exp(ygrid); % Convert log process back to levels
            par.pmat = pmat; % Transition matrix.
        end
                
        %% Tauchen's Method
        function [y, pi] = tauchen(mu, rho, sigma, N, m)
            %% Construct equally spaced grid.
            ar_mean = mu / (1 - rho); % Mean of stationary AR(1) process
            ar_sd = sigma / sqrt(1 - rho^2); % Std. dev of stationary AR(1) process
            
            y1 = ar_mean - (m * ar_sd); % Smallest grid point
            yn = ar_mean + (m * ar_sd); % Largest grid point
            
            y = linspace(y1, yn, N); % Equally spaced grid.
            d = y(2) - y(1); % Step size.
	        
            %% Compute transition probability matrix
            ymatk = repmat(y, N, 1); % States next period
            ymatj = mu + rho * ymatk'; % States this period
        
	        pi = normcdf(ymatk, ymatj - (d/2), sigma) - normcdf(ymatk, ymatj + (d/2), sigma);
	        pi(:, 1) = normcdf(y(1), mu + rho*y - (d/2), sigma);
	        pi(:, N) = 1 - normcdf(y(N), mu + rho*y + (d/2), sigma);
        end
        
        %% Utility function.
        function u = utility(c, n, par)

            %% Leisure.
            un = ((1-n).^(1+(1/par.nu)))./(1+(1/par.nu)); 

            %% CRRA utility.
            if par.sigma == 1
                
                uc = log(c); % Log utility of consumption.
            else
                uc = (c.^(1-par.sigma))./(1-par.sigma); % CRRA utility.
            end
            
             % Total.
            u = uc + par.gamma*un;

        end
    end
end
