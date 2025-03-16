%% File Info.
%{
    solve.m
    -------
    This code solves the model using Value Function Iteration (VFI).
%}

%% Solve class.

classdef solve
    methods(Static)
        
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
            g0 = tauk * r * kmat + taun * w * n - delta * tauk * kmat; % Government spending.
            c0 = max(y0 - i0 - g0, 1e-6); % Ensure non-negative consumption.
            v0 = model.utility(c0, g0, par) / (1 - beta); % Initial guess.

            v1 = zeros(klen, Alen); % Value function.
            k1 = zeros(klen, Alen); % Policy function for k'.
            g1 = zeros(klen, Alen); % Government spending policy function.

            crit = 1e-6;
            maxiter = 10000;
            diff = 1;
            iter = 0;
            
            fprintf('------------Starting Value Function Iteration------------\n\n')
            
            while diff > crit && iter < maxiter  
                
                for j = 1:Alen  
                    for p = 1:klen  

                        % Compute output.
                        y = Agrid(j) * (kgrid(p)^alpha);
                        
                        % Solve maximization problem.
                        ev = v0 * pmat(j,:)';  
                        vall = -inf(klen, 1); % Initialize value function array

                        for ind = 1:klen
                            % Compute investment choices dynamically.
                            i = max(kgrid(ind) - (1 - delta) * kgrid(p), 1e-6);

                            % Compute government spending dynamically.
                            g = tauk * r * kgrid(ind) + taun * w * n - delta * tauk * kgrid(p);

                            % Compute consumption from budget constraint.
                            c = max(y - i - g, 1e-6); % Ensure positive consumption.

                            % Compute value function.
                            vall(ind) = model.utility(c, g, par) + beta * ev(ind);
                        end
                        
                        % Choose the best policy.
                        [vmax, ind] = max(vall);
                        

                        % Store values.
                        v1(p,j) = vmax;
                        k1(p,j) = kgrid(min(max(ind, 1), klen)); % Keep within bounds
                        g1(p,j) = g;
                    end
                end
                
                % Check for convergence.
                diff = norm(v1 - v0);
                v0 = v1;
                iter = iter + 1;
                
                if mod(iter, 25) == 0
                    fprintf('Iteration: %d.\n', iter);
                end
            end
                
            fprintf('\nConverged in %d iterations.\n\n', iter)
            fprintf('------------End of Value Function Iteration------------\n')
            
            %% Store results.
            sol.y = Amat .* (kmat.^alpha);
            sol.k = k1;
            sol.i = k1 - (1 - delta) * kmat; 
            sol.g = tauk * r * k1+ taun * w * n - delta * tauk * k1;
            sol.c = sol.y - sol.i - sol.g; 
            sol.v = v1;

            
            
        end
    end
end