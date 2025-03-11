%% File Info.

%{
    simulate.m
    ----------
    This code simulates the life cycle model with the optimal policies.
%}

%% Simulate class.

classdef simulate
    methods(Static)
        %% Simulate the model. 
        
        function sim = grow(par, sol)            
            %% Set up.
            
            agrid = par.agrid;    % Asset grid (state variable).
            T = par.T;            % Time periods.
            tr = par.tr;          % Retirement period
            r = par.r;            % Interest rate
            kappa = par.kappa;    % Pension fraction
            ybar = par.ybar;      % Fixed income in working years
            
            a_pol = sol.a;        % Policy function for assets.
            c_pol = sol.c;        % Policy function for consumption.
            
            a = zeros(T, 1);      % Container for asset holdings
            c = zeros(T, 1);      % Container for consumption
            y = zeros(T, 1);      % Container for income
            u = zeros(T, 1);      % Container for utility
            
            %% Begin simulation from period 1 with initial wealth.
            
            a(1) = par.a0;        % Initial asset level from parameters
            
            %% Simulate entire lifecycle
            
            for t = 1:T
                % Determine income based on working/retirement period
                if t < tr
                    y(t) = ybar;  % Working period income
                else
                    y(t) = kappa * ybar;  % Retirement income (pension)
                end
                
                % Find index of current asset level in the grid
                [~, a_idx] = min(abs(a(t) - agrid));
                
                % Get consumption for current period
                c(t) = c_pol(a_idx, t);
                
                % Calculate utility
                u(t) = model.utility(c(t), par);
                
                % Update assets for next period if not the last period
                if t < T
                    a(t+1) = a_pol(a_idx, t);
                end
            end
            
            %% Store results in sim structure
            
            sim = struct();
            sim.a = a;      % Asset holdings at each age
            sim.c = c;      % Consumption at each age
            sim.y = y;      % Income at each age
            sim.u = u;      % Utility at each age
        end
    end
end