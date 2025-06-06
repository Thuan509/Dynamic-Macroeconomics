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
        
        function sol = grow(par)            
            %% Structure array for model solution.
            
            sol = struct();
            
            %% Model parameters, grids, and functions.
            beta = par.beta; % Discount factor.
            alpha = par.alpha; % Capital share.
            delta = par.delta; % Depreciation.
            tauk = par.tauk; % Capital tax rate.
            taun = par.taun; % Labor tax rate.
            w = par.w; % Wage.
            r = par.r; % Return on capital.
            n = par.n; % Labor supply

            klen = par.klen; % Grid size for k.
            kgrid = par.kgrid; % Grid for k.
            Alen = par.Alen; % Grid size for A.
            Agrid = par.Agrid; % Grid for A.
            pmat = par.pmat; % Transition probabilities for A.

            kmat = repmat(kgrid,1,Alen); % k for each value of A.
            Amat = repmat(Agrid,klen,1); % A for each value of k.

            
            %% Initialize value function.
            y0 = Amat .* (kmat.^alpha); % Output.
            i0 = kmat - (1 - delta) * kmat; % Investment.
            g0 = tauk * r * kmat + taun * w * n - delta * tauk * kmat; % Initial government spending estimate.
            c0 = max(y0 - i0 - g0, 1e-6); % Ensure non-negative consumption.
            v0 = model.utility(c0, g0, par) / (1 - beta); % Initial guess.

            v1 = zeros(klen, Alen); % Value function.
            k1 = zeros(klen, Alen); % Policy function for k'.
            g1 = zeros(klen, Alen); % Government spending policy function.
                        
            crit = 1e-6;
            maxiter = 10000;
            diff = 1;
            iter = 0;
            
            fprintf('------------Beginning Value Function Iteration.------------\n\n')
            
            while diff > crit && iter < maxiter % Iterate on the Bellman Equation until convergence.
                
                for p = 1:klen % Loop over the k-states.
                    for j = 1:Alen % Loop over the A-states.

                        % Compute output.
                        y = Agrid(j) * (kgrid(p)^alpha);

                        % Solve maximization problem over investment and government spending.
                        ev = v0 * pmat(j,:)';  % Expected value function over future A states.
                        vall = -inf(klen, 1); % Initialize value function array
                        g_opt = zeros(klen, 1); % Store optimal g values

                        for ind = 1:klen
                            % Compute investment choices dynamically.
                            i = kgrid(ind) - (1 - delta) * kgrid(p);

                            % Compute government spending dynamically.
                            g = max(tauk * r * kgrid(ind) + taun * w * n - delta * tauk * kgrid(p), 1e-6);
                            g_opt(ind) = g;

                            % Compute consumption from budget constraint.
                            c = max(y - i - g, 1e-6);

                            % Compute value function.
                            vall(ind) = model.utility(c, g, par) + beta * ev(ind);
                        end
                        
                        % Choose the best policy.
                        [vmax, ind] = max(vall);
                        
                        % Store values.
                        v1(p,j) = vmax;
                        k1(p,j) = kgrid(min(max(ind, 1), klen)); % Keep within bounds
                        g1(p,j) = g_opt(ind); % Store optimal government spending
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
            
            %% Store results.
            
            sol.y = Amat.*(kmat.^alpha);                                         % Output.
            sol.k = k1;                                                                       % Capital policy function.
            sol.i = k1 - ((1 - delta).* kmat);                                       % Investment policy function.
            sol.g = g1; % Optimized government spending policy function.
            sol.c = sol.y - sol.i - sol.g;                                               % Consumption policy function.
            sol.v = v1;                                                                       % Value function.
        end
    end
end
