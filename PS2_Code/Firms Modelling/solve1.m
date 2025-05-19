%% File Info.

%{

    solve1.m
    -------
    This code solve the model.

%}

%% Solve class.

classdef solve1
    methods(Static)
        %% Solve the model using VFI.
        
        function sol = firm_problem(par)
            %% Structure array for model solution.
            sol = struct();

            %% Model parameters and grids
            beta = par.beta;
            delta = par.delta;
            wt = par.wt;
            alpha = par.alpha;

            klen = par.klen;
            kgrid = par.kgrid;
            Alen = par.Alen;
            Agrid = par.Agrid;
            pmat = par.pmat;
            plen = par.plen;
            pgrid = par.pgrid;
            mtrix = par.mtrix;

            %% Initialization
            v0 = zeros(klen,Alen,plen);
            v1 = nan(klen,Alen,plen);
            k1 = nan(klen,Alen,plen);
            r1 = nan(klen,Alen,plen);
            e1 = nan(klen,Alen,plen);
            p1 = nan(klen,Alen,plen);
            x1 = nan(klen,Alen,plen);

            crit = 1e-6;
            maxiter = 10000;
            diff = 1;
            iter = 0;

            fprintf('------------Beginning Value Function Iteration.------------\n\n')

            while diff > crit && iter < maxiter

                for p = 1:klen % Loop over capital states
                    for j = 1:Alen % Loop over productivity states
                        for i = 1:plen % Loop over price states

                            % Macro variables
                            xt = (((1-alpha) * Agrid(j) * (kgrid(p)^alpha)) / wt)^(1/alpha); % Optimal input choice
                            rev = model1.production(Agrid(j), kgrid(p), xt, par); % Revenue
                            exp = model1.total_cost(xt, kgrid(p), pgrid(i), par); % Total cost
                            prof = rev - exp; % Profit

                            % Expected continuation value
                            ev = zeros(klen,1);
                            for q = 1:plen
                                ev = ev + squeeze(v0(:,:,q)) * (pmat(j,:)') * mtrix(i,q); 
                            end

                            % Bellman choice
                            vall = prof + beta * ev; % profit + discounted expected future value

                            [vmax, ind] = max(vall); % Find best next period capital

                            % Store optimal policies
                            v1(p,j,i) = vmax;
                            k1(p,j,i) = kgrid(ind); % Optimal k'
                            r1(p,j,i) = rev;
                            e1(p,j,i) = exp(ind);
                            p1(p,j,i) = prof(ind);
                            x1(p,j,i) = xt;
                        end
                    end
                end

                diff = max(abs(v1(:) - v0(:))); % convergence check
                v0 = v1;
                iter = iter + 1;

                if mod(iter, 25) == 0
                    fprintf('Iteration: %d\n', iter);
                end
            end

            fprintf('\nConverged in %d iterations.\n\n', iter)
            fprintf('------------End of Value Function Iteration.------------\n')

            %% Save output
            sol.v = v1;
            sol.k = k1;
            sol.r = r1;
            sol.e = e1;
            sol.p = p1;
            sol.x = x1;

        end
    end
end
