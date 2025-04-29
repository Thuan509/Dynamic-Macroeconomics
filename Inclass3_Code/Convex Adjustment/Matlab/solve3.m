%% File Info.
%{
    solve3.m
    -------
    This code solves the model with discrete investment choice and borrowing constraint.
%}

%% Solve class.

classdef solve3
    methods(Static)
        %% Solve the model using VFI. 
        
        function sol = firm_problem(par)            
            %% Structure array for model solution.
            
            sol = struct();
            
            %% Model parameters, grids, and functions.
            
            beta = par.beta; % Discount factor.
            delta = par.delta; % Depreciation rate.
            kappa = par.kappa; % Fixed cost parameter (κ)

            klen = par.klen; % Grid size for k.
            kgrid = par.kgrid; % Grid for k (state and choice).

            Alen = par.Alen; % Grid size for A.
            Agrid = par.Agrid; % Grid for A.
            pmat = par.pmat; % Transition matrix for A.
            lambda = par.lambda; % Proportion of revenue kept when investing

            %% Value Function Iteration.

            v0 = zeros(klen,Alen); % Initial guess of value function.

            v1 = nan(klen,Alen); % Container for value function.
            k1 = nan(klen,Alen); % Container for K'.
            i1 = nan(klen,Alen); % Container for investment.
            r1 = nan(klen,Alen); % Container for revenue.
            e1 = nan(klen,Alen); % Container for investment expenditure.
            p1 = nan(klen,Alen); % Container for profit.
            invest_decision = nan(klen,Alen); % Store discrete choice (1=Invest, 0=No Invest).

            crit = 1e-6;
            maxiter = 10000;
            diff = 1;
            iter = 0;
            
            fprintf('------------ Beginning Value Function Iteration ------------\n\n')
            
            while diff > crit && iter < maxiter % Iterate on the Bellman Equation until convergence.
                
                for p = 1:klen % Loop over the K-states.
                    for j = 1:Alen % Loop over the A-states.

                        % Compute revenue
                        rev = model2.production(Agrid(j), kgrid(p), par); 
                        r1(p,j) = rev;  % Store revenue

                        % Case 1: No Investment
                        k_no_invest = (1 - delta) * kgrid(p); % Capital next period if no investment.
                        profit_no_invest = rev; % Full revenue, no investment cost.
                        ev_no_invest = 0;
                        % Find the index of k_no_invest in kgrid or the closest value
                        [~, idx] = min(abs(kgrid - k_no_invest));
                        ev_no_invest = v0(idx,:) * pmat(j,:)'; % Expected future value.
                        value_no_invest = profit_no_invest + beta * ev_no_invest; 

                        % Case 2: Invest
                        value_invest = -Inf; % Default to -Inf in case no profitable investment exists.
                        best_k_prime = k_no_invest; % Default next-period capital.
                        best_exp_invest = 0;
                        best_invest = 0;

                        for q = 1:klen  % Loop over possible K' choices
                            k_prime = kgrid(q);
                            investment = k_prime - (1 - delta) * kgrid(p); % Compute investment amount
                            
                            % Skip if investment exceeds borrowing constraint (it ≤ AKα)
                            if investment > rev
                                continue;
                            end
                            
                            % Calculate non-convex adjustment cost: C(K',A,K) = κK
                            cost = kappa * kgrid(p);
                            
                            % Total investment cost including the price of investment
                            if par.p > 0 % If price of investment is specified
                                exp_invest = cost + par.p * investment;
                            else
                                exp_invest = cost;
                            end
                            
                            % Revenue after loss and profit
                            revenue_after_loss = lambda * rev;
                            profit_invest = revenue_after_loss - exp_invest;
                            
                            if profit_invest >= 0 % Only consider cases where investment is profitable
                                ev_invest = v0(q,:) * pmat(j,:)'; % Expected future value with new capital
                                temp_value = profit_invest + beta * ev_invest; % Compute total value.

                                % Keep the best investment option
                                if temp_value > value_invest
                                    value_invest = temp_value;
                                    best_k_prime = k_prime;
                                    best_exp_invest = exp_invest;
                                    best_invest = investment;
                                end
                            end
                        end

                        % Choose between Investing and Not Investing
                        if value_invest > value_no_invest
                            v1(p,j) = value_invest; % Choose to invest
                            k1(p,j) = best_k_prime;
                            i1(p,j) = best_invest;
                            e1(p,j) = best_exp_invest;
                            p1(p,j) = lambda * rev - best_exp_invest; % Profit after investment with revenue loss
                            invest_decision(p,j) = 1; % Firm chooses to invest
                        else
                            v1(p,j) = value_no_invest; % Choose not to invest
                            k1(p,j) = k_no_invest;
                            i1(p,j) = 0; % No investment.
                            e1(p,j) = 0; % No expenditure.
                            p1(p,j) = profit_no_invest; % Keep full revenue.
                            invest_decision(p,j) = 0; % Firm does not invest
                        end

                    end
                end
                
                diff = norm(v1 - v0); % Check for convergence.
                v0 = v1; % Update value function guess.
                
                iter = iter + 1; % Update counter.
                
                % Print progress
                if mod(iter, 25) == 0
                    fprintf('Iteration: %d, Diff: %f\n', iter, diff)
                end

            end
                
            fprintf('\nConverged in %d iterations with diff = %f.\n\n', iter, diff)
            
            fprintf('------------ End of Value Function Iteration ------------\n')
            
            %% Store results in solution structure.
            
            sol.v = v1; % Value function.
            sol.k = k1; % Optimal capital policy.
            sol.i = i1; % Investment decision.
            sol.r = r1; % Revenue function.
            sol.e = e1; % Investment cost function.
            sol.p = p1; % Profit function.
            sol.invest_decision = invest_decision; % Discrete choice indicator.
            
        end
        
    end
end