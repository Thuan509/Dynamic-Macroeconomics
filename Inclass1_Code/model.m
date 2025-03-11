%% File Info.

%{

    model.m
    -------
    This code sets up the Life Cycle Model parameters.

%}

%% Model class.

classdef model
    methods(Static)
        %% Set up structure array for model parameters and set the simulation parameters.
        
        function par = setup()            
            %% Structure array for model parameters.
            
            par = struct();
            
            %% Preferences.
            
            par.beta = 0.96; % Discount factor: Lower values of this mean that consumers are impatient and consume more today.
            par.sigma = 2.00; % CRRA: Higher values of this mean that consumers are risk averse and do not want to consume too much today.
            
            assert(par.beta > 0 && par.beta < 1.00,'Discount factor should be between 0 and 1.\n')
            assert(par.sigma > 0,'CRRA should be at least 0.\n')

            %% Income and Technology Parameters

            par.ybar = 10; % Fixed labor income during working years
            par.kappa = 0.3; % Pension fraction of income (0 < kappa < 1)
            par.r = 0.1; % Fixed interest rate
            par.a0 = 5; % Initial wealth endowment (a0 > 0)
            par.aT = 0; % Terminal condition (no bequests)

            assert(par.ybar > 0,'The fixed income should be larger than 0.\n')
            assert(par.kappa > 0 && par.kappa < 1,'The pension fraction should be between 0 and 1.\n')
            assert(par.a0 > 0,'The initial endowment should be greater than 0.\n')
            assert(par.aT == 0,'The terminal wealth should be equal to 0.\n')

            %% Simulation parameters.

            par.seed = 2025; % Seed for simulation.
            par.T = 4; % Total lifespan (4 periods)
            par.tr = 2; % Retirement begins at period 2 (work periods: 0,1; retirement periods: 2,3)
        end
        
        %% Generate state grids.
        
        function par = gen_grids(par)
            %% Wealth grid.
            par.alen = 300; % Grid size for a.
            par.amax = 30; % Upper bound for a.
            par.amin = 0; % Minimum a.
            
            assert(par.alen > 5,'Grid size for a should be positive and greater than 5.\n')
            assert(par.amax > par.amin,'Minimum a should be less than maximum value.\n')
            
            par.agrid = linspace(par.amin,par.amax,par.alen)'; % Equally spaced, linear grid for a.
        end
        
        %% Utility function.
        
        function u = utility(c,par)
            %% CRRA utility.
            
            if par.sigma == 1
                u = log(c); % Log utility.
            else
                u = (c.^(1-par.sigma))./(1-par.sigma); % CRRA utility.
            end
        end
        
    end
end