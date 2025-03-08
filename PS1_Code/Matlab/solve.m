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
            sigma = par.sigma; % CRRA.
            tauk = par.tauk; % Capital tax rate.
            taun = par.taun; % Labor tax rate.
            g = par.g;  % Government spending.
            t = par.t; % Transfers.
            w = par.w; % Wage.
            r = par.r; % Return on capital.

            klen = par.klen; % Grid size for k.
            kgrid = par.kgrid; % Grid for k.
            Alen = par.Alen; % Grid size for A.
            Agrid = par.Agrid; % Grid for A.
            pmat = par.pmat; % Transition probabilities for A.

            kmat = repmat(kgrid,1,Alen); % k for each value of A.
            Amat = repmat(Agrid,klen,1); % A for each value of k.
            
            %% Initialize value function.
            y0 = Amat .* (kmat.^alpha); % Output.
            i0 = (1 - delta) * kmat; % Depreciated capital.
            c0 = y0 - i0; % Consumption.
            v0 = model.utility(c0, g, par) / (1 - beta); % Initial guess.

            v1 = zeros(klen, Alen); % Value function.
            k1 = zeros(klen, Alen); % Policy function for k'.

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
                        
                        % Compute investment choices.
                        i = kgrid - (1 - delta) * kgrid(p);
                        
                        % Compute consumption from budget constraint.
                        total_income = (1 - tauk) * r * kgrid + (1 - taun) * w + delta * tauk * kgrid + t;
                        c = total_income - i;

                        % Solve maximization problem.
                        ev = v0 * pmat(j,:)';  
                        vall = model.utility(c, g, par) + beta * ev;  
                        vall(c < 0) = -inf;  
                        [vmax, ind] = max(vall);  

                        % Store values.
                        v1(p,j) = vmax;
                        k1(p,j) = kgrid(ind);
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
            sol.c = sol.y - sol.i;
            sol.v = v1;
            
        end
    end
end