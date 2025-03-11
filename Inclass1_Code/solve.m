%% File Info.

%{
    solve.m
    -------
    This code solves the Life Cycle Model using backward iteration.
%}

%% Solve class.

classdef solve
    methods(Static)
        %% Solve the life cycle model using backward iteration.
        
        function sol = grow(par)            
            %% Structure array for model solution.
            
            sol = struct();

            %% Model parameters, grids, and functions.
            
            beta = par.beta; % Discount factor
            r = par.r; % Interest rate
            kappa = par.kappa; % Pension fraction
            ybar = par.ybar; % Fixed income in working years
            tr = par.tr; % Retirement period
            T = par.T; % Total lifespan

            alen = par.alen; % Grid size for assets
            agrid = par.agrid; % Grid for assets
            
            %% Initialize Value and Policy Functions
            
            v = zeros(alen, T); % Value function
            a_pol = zeros(alen, T); % Policy function for next period assets (a')
            c_pol = zeros(alen, T); % Policy function for consumption

            fprintf('------------ Solving Life Cycle Model Backward ------------\n\n')

            %% Solve backward from last period to first
            
            for t = T:-1:1  % Start from last period (T) and move backward to first period (1)
                fprintf('Solving for period %d...\n', t)
                
                % Determine income based on working/retirement period
                if t < tr
                    y_t = ybar; % Working period: regular labor income
                else
                    y_t = kappa * ybar; % Retirement period: pension income
                end
                
                for i_a = 1:alen  % Loop over current asset states
                    a_t = agrid(i_a);  % Current assets
                    
                    if t == T  % Last period: consume everything (aT = 0)
                        % In the last period, optimal to consume all wealth
                        a_pol(i_a, t) = par.aT;  % Terminal condition: aT = 0
                        c_pol(i_a, t) = (1 + r) * (a_t + y_t);  % Consume everything
                        v(i_a, t) = model.utility(c_pol(i_a, t), par);  % Terminal value
                    else
                        % For periods before the last one, find optimal a' (savings)
                        value_vec = zeros(alen, 1);  % Store values for each possible a'
                        
                        for i_ap = 1:alen  % Loop over all possible next period assets
                            a_next = agrid(i_ap);  % Possible next period assets
                            
                            % Consumption based on budget constraint: c = (1+r)(a+y) - a'
                            c_t = (1 + r) * (a_t + y_t) - a_next;
                            
                            if c_t > 0  % Only consider positive consumption
                                % Value = utility today + discounted future value
                                value_vec(i_ap) = model.utility(c_t, par) + beta * v(i_ap, t+1);
                            else
                                value_vec(i_ap) = -inf;  % Negative consumption not allowed
                            end
                        end
                        
                        % Find the optimal a' that maximizes value
                        [v_max, i_a_max] = max(value_vec);
                        a_next_opt = agrid(i_a_max);
                        
                        % Store optimal values and policies
                        v(i_a, t) = v_max;
                        a_pol(i_a, t) = a_next_opt;
                        c_pol(i_a, t) = (1 + r) * (a_t + y_t) - a_next_opt;
                    end
                end
            end

            fprintf('------------ Solution Completed ------------\n')

            %% Store Solution
            
            sol.a = a_pol;  % Asset policy function
            sol.v = v;      % Value function
            sol.c = c_pol;  % Consumption policy function
        end
    end
end