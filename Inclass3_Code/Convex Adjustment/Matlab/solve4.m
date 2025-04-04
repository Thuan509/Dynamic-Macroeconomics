%% File Info.

%{

    solve.m
    -------
    This code solves the model.

%}

%% Solve class.

classdef solve
    methods(Static)
        %% Solve the model using VFI. 
        
        function sol = firm_problem(par)            
            %% Structure array for model solution.
            
            sol = struct();
            
            %% Model parameters, grids and functions.
            
            beta = par.beta; % Discount factor.
            delta = par.delta; % Depreciation rate.

            klen = par.klen; % Grid size for k.
            kgrid = par.kgrid; % Grid for k (state and choice).

            Alen = par.Alen; % Grid size for A.
            Agrid = par.Agrid; % Grid for A.
            pmat = par.pmat; % Transition matrix for A.

            %% Value Function iteration.

            v0 = zeros(klen,Alen); % Guess of value function is zero profit.

            v1 = nan(klen,Alen); % Container for V.
            k1 = nan(klen,Alen); % Container for K'.
            i1 = nan(klen,Alen); % Container for i.
            r1 = nan(klen,Alen); % Container for revenue.
            e1 = nan(klen,Alen); % Container for investment expenditure.
            p1 = nan(klen,Alen); % Container for profit.

            crit = 1e-6;
            maxiter = 10000;
            diff = 1;
            iter = 0;
            
            fprintf('------------Beginning Value Function Iteration.------------\n\n')
            
            while diff > crit && iter < maxiter % Iterate on the Bellman Equation until convergence.
                
                for p = 1:klen % Loop over the K-states.
                    for j = 1:Alen % Loop over the A-states.

                        % Macro variables.
                        rev = model.production(Agrid(j),kgrid(p),par); % Revenue given A and K.
                        exp = model.total_cost(kgrid(p),par); % Total investment expenditure given K.
                        prof = rev-exp; % Profit.
                        invest = par.kgrid-(1-delta).*kgrid(p); % Investment in new capital.
                         % Apply borrowing constraint: i_t \leq A K^alpha
                        %borrow_limit = Agrid(j) * (kgrid(p) ^ alpha);
                        %invest(invest > borrow_limit) = -inf; % Restrict investment choices.

                        % Solve the maximization problem.
                        ev = v0*pmat(j,:)'; %  The next-period value function is the expected value function over each possible next-period A, conditional on the current state j.
                        vall = prof + beta*ev; % Compute the value function for each choice of K', given K.
                        vall(prof<0) = -inf; % Set the value function to negative infinity when profit < 0.
                        vall(invest> rev) = -inf; % Set the value function to negative infinity when profit < 0.
                         vall(invest<0) = -inf; 
                        [vmax,ind] = max(vall); % Maximize: vmax is the maximized firm value; ind is where it is in the grid.
                    
                        % Store values.
                        v1(p,j) = vmax; % Maximized firm value.
                        k1(p,j) = kgrid(ind); % Optimal K'.
                        i1(p,j) = invest(ind); % Optimal i.
                        r1(p,j) = rev; % Total revenue.
                        e1(p,j) = exp(ind); % Total cost.
                        p1(p,j) = prof(ind); % Profits.

                    end
                end
                
                diff = norm(v1-v0); % Check for convergence.
                v0 = v1; % Update guess of v.
                
                iter = iter + 1; % Update counter.
                
                % Print counter.
                if mod(iter,25) == 0
                    fprintf('Iteration: %d.\n',iter)
                end

            end
                
            fprintf('\nConverged in %d iterations.\n\n',iter)
            
            fprintf('------------End of Value Function Iteration.------------\n')
            
            %% Macro variables, value, and policy functions.
            
            sol.v = v1; % Firm value.
            sol.k = k1; % Capital policy function.
            sol.i = i1; % Investment policy function.
            sol.r = r1; % Revenue function.
            sol.e = e1; % Investment expenditure function.
            sol.p = p1; % Profit function.
            
        end
        
    end
end