%% File Info.

%{
    solve.m
    -------
    This code solves the Life Cycle Model with stochastic income using backward iteration.
%}

%% Solve class.

classdef solve2
    methods(Static)
        %% Solve the life cycle model with stochastic income using backward iteration.
        
        function sol = grow(par)            
            %% Structure array for model solution.
            
            sol = struct();

            %% Model parameters
            
            beta = par.beta;      % Discount factor
            r = par.r;            % Interest rate
            kappa = par.kappa;    % Pension fraction
            tr = par.tr;          % Retirement age
            T = par.T;            % Total lifespan

            alen = par.alen;      % Number of asset grid points
            agrid = par.agrid;    % Asset grid
            ylen = par.ylen;      % Number of income states
            ygrid = par.ygrid;    % Income grid
            ytrans = par.ytrans;  % Income transition matrix
            
            %% Initialize Value and Policy Functions
            
            % Value function: v(a, y, t)
            v = zeros(alen, ylen, T);
            
            % Policy functions: a'(a, y, t) and c(a, y, t)
            a_pol = zeros(alen, ylen, T);
            c_pol = zeros(alen, ylen, T);

            fprintf('------------ Solving Stochastic Life Cycle Model Backward ------------\n\n')

            %% Last period solution (consume everything)
            
            % In the last period, optimal policy is to consume all resources
            for i_y = 1:ylen
                if T >= tr
                    % Retirement income (pension based on last working period income)
                    % For simplicity, we use average income as the last working period income
                    y_T = kappa * mean(ygrid);
                else
                    % Still working in last period (should not happen with T > tr)
                    y_T = ygrid(i_y);
                end
                
                for i_a = 1:alen
                    a_T = agrid(i_a);
                    
                    % Terminal condition: consume everything
                    c_T = (1 + r) * (a_T + y_T);
                    
                    % Store policy and value
                    a_pol(i_a, i_y, T) = par.aT;  % No savings in last period
                    c_pol(i_a, i_y, T) = c_T;
                    v(i_a, i_y, T) = model.utility(c_T, par);
                end
            end
            
            % Last working period income (needed for pension calculation)
            last_work_income = zeros(alen, ylen);
            
            %% Solve backward from period T-1 to period 1
            
            for t = (T-1):-1:1  % Work backward
                fprintf('Solving for period %d...\n', t)
                
                % Check if this is retirement or working period
                is_retired = (t >= tr);
                
                % This will track the last working period income for retirement
                if t == tr-1
                    is_last_work = true;
                else
                    is_last_work = false;
                end
                
                for i_y = 1:ylen  % Loop over income states
                    y_t = ygrid(i_y);  % Current income
                    
                    % In retirement, income depends on last working period income
                    if is_retired
                        % For each asset level, income is constant in retirement
                        % (based on last working period)
                        if t == tr
                            % First retirement period: pension depends on current y
                            income = kappa * y_t;
                        else
                            % Later retirement periods: pension is constant
                            income = kappa * y_t;  % Simplified
                        end
                    else
                        % Working period: income is stochastic
                        income = y_t;
                    end
                    
                    for i_a = 1:alen  % Loop over asset states
                        a_t = agrid(i_a);  % Current assets
                        
                        % Total resources available
                        resources = (1 + r) * (a_t + income);
                        
                        % Optimization for each current state (a_t, y_t)
                        value_vec = zeros(alen, 1);  % Value for each possible a'
                        
                        for i_ap = 1:alen  % Loop over possible next-period assets
                            a_next = agrid(i_ap);  % Next period assets
                            
                            % Consumption from budget constraint
                            c_t = resources - a_next;
                            
                            if c_t <= 0
                                % Negative consumption not allowed
                                value_vec(i_ap) = -inf;
                            else
                                % Calculate expected continuation value
                                expected_value = 0;
                                
                                if is_retired
                                    % In retirement, income is deterministic
                                    expected_value = v(i_ap, i_y, t+1);
                                else
                                    % In working periods, income follows Markov process
                                    for i_yp = 1:ylen
                                        expected_value = expected_value + ytrans(i_y, i_yp) * v(i_ap, i_yp, t+1);
                                    end
                                end
                                
                                % Value = current utility + discounted future value
                                value_vec(i_ap) = model.utility(c_t, par) + beta * expected_value;
                            end
                        end
                        
                        % Find optimal a' and associated value
                        [v_max, i_ap_max] = max(value_vec);
                        a_next_opt = agrid(i_ap_max);
                        c_opt = resources - a_next_opt;
                        
                        % Store optimal policy and value
                        v(i_a, i_y, t) = v_max;
                        a_pol(i_a, i_y, t) = a_next_opt;
                        c_pol(i_a, i_y, t) = c_opt;
                        
                        % Store last working period income (for pension calculation)
                        if is_last_work
                            last_work_income(i_a, i_y) = income;
                        end
                    end
                end
            end

            fprintf('------------ Solution Completed ------------\n')

            %% Store Solution
            
            sol.a = a_pol;        % Asset policy function a'(a, y, t)
            sol.c = c_pol;        % Consumption policy function c(a, y, t)
            sol.v = v;            % Value function v(a, y, t)
            sol.last_work_income = last_work_income;  % Income in last working period
        end
    end
end