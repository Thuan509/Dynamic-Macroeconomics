%% File Info.

%{

    solve1.m
    -------
    This code solves the model.

%}

%% Solve class.

classdef solve1
    methods(Static)
        %% Solve the model using VFI. 
        
        function sol = firm_problem(par)            
            %% Structure array for model solution.
            
            sol = struct();
            
            %% Model parameters, grids and functions.
            
            beta = par.beta; % Discount factor.
            delta = par.delta; % Depreciation rate.

            wt = par.wt;  % Variable costs.
            %gamma = par.gamma; 
            alpha = par.alpha;
            %sigma = par.sigma;

            klen = par.klen; % Grid size for k.
            kgrid = par.kgrid; % Grid for k (state and choice).

            Alen = par.Alen; % Grid size for A.
            Agrid = par.Agrid; % Grid for A.
            pmat = par.pmat; % Transition matrix for A.

            plen = par.plen; % Grid size for p.
            pgrid = par.pgrid; % Grid for p.
            mtrix = par.mtrix; % Transition matrix for p.

            %% Value Function iteration.

            v0 = zeros(klen,Alen, plen); % Guess of value function is zero profit.

            v1 = nan(klen,Alen, plen); % Container for V.
            k1 = nan(klen,Alen, plen); % Container for K'.
            i1 = nan(klen,Alen, plen); % Container for i.
            r1 = nan(klen,Alen, plen); % Container for revenue.
            e1 = nan(klen,Alen, plen); % Container for investment expenditure.
            p1 = nan(klen,Alen, plen); % Container for profit.
            x1 = nan(klen,Alen, plen); % Container input cost. 

            crit = 1e-6;
            maxiter = 10000;
            diff = 1;
            iter = 0;
            
            fprintf('------------Beginning Value Function Iteration.------------\n\n')
            
            while diff > crit && iter < maxiter % Iterate on the Bellman Equation until convergence.
                
                for p = 1:klen % Loop over the K-states.
                    for j = 1:Alen % Loop over the A-states.
                        for i = 1:plen % Loop over p-states.

                            % Macro variables.
                            xt = (((1-alpha).*Agrid(j).*(kgrid(p)^alpha))/wt)^(1./alpha); % Define the optimal input choice.
                            rev = model1.production(Agrid(j),kgrid(p),xt, par); % Revenue given A and K.
                            exp = model1.total_cost(xt,kgrid(p),pgrid(i),par); % Total investment expenditure given K.
                            prof = rev-exp; % Profit.
                            invest = pgrid(i).* (par.kgrid-(1-delta).*kgrid(p)); % Investment in new capital.

                            ev = zeros(klen,1);

                            for q =1:plen    
                                % Solve the maximization problem.
                                ev = ev + squeeze(v0(:,:,q)) * (pmat(j,:)') * mtrix(i,q); %  The next-period value function is the expected value function over each possible next-period A, conditional on the current state j.
                            end

                           vall = prof + beta*ev; % Compute the value function for each choice of K', given K.
                           vall(prof<0) = -inf; % Set the value function to negative infinity when profit < 0.
                            vall(invest> rev) = -inf; % Set the value function to negative infinity when profit < 0.
                            vall(invest<0) = -inf; 

                            
                            [vmax,ind] = max(vall); % Maximize: vmax is the maximized firm value; ind is where it is in the grid.
                        
                            % Store values.
                            v1(p,j,i) = vmax; % Maximized firm value.
                            k1(p,j,i) = kgrid(ind); % Optimal K'.
                            i1(p,j,i) = invest(ind); % Optimal i.
                            r1(p,j,i) = rev; % Total revenue.
                            e1(p,j,i) = exp(ind); % Total cost.
                            p1(p,j,i) = prof(ind); % Profits.
                            x1(p,j,i) = xt(ind); % Input costs.

                        end
                    end
                end
                
                %diff = norm(v1-v0); % Check for convergence.
                diff = max(abs(v1(:)-v0(:))); % Check for convergence.
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
            sol.x = x1; % Input cost function.
            
        end
        
    end
end