%% File Info.

%{

    solve.m
    -------
    % This code solves the life cycle model by computing optimal consumption and asset choices.

%}

%% Solve class.

classdef solve
    methods(Static)
        function sol = lc(par)
            sol = struct();
            T = par.T; % Number of periods
            tr = par.tr; % Retirement period
            beta = par.beta; % Discount factor
            alen = par.alen; % Asset grid size
            agrid = par.agrid; % Asset grid
            ylen = par.ylen; % Income states
            ygrid = par.ygrid; % Income grid
            pmat = par.pmat; % Transition probability matrix
            r_values = [-0.5, 0.5]; % Range of risk-free rates
            kappa = par.kappa; % Retirement income multiplier

            alpha=linspace(0.0,1.0,100); %Grid space for alpha

            r=par.r;
            % Set borrowing limit (default if missing)
            if ~isfield(par, 'borrow_limit')
                par.borrow_limit = 0.5; % Default limit
            end
            
            amat = repmat(agrid,1,100);

            v1 = nan(alen, T, ylen, 2); % Value function
            a1 = nan(alen, T, ylen, 2); % Asset policy function
            c1 = nan(alen, T, ylen, 2); % Consumption policy function
            alpha1 = nan(alen, T, ylen, 2); %alpha policy function
            
            for age = 1:T
                for rr = 1:2
                    if T - age + 1 == T
                        c1(:, T, :,:) = repmat(reshape(agrid + kappa * ygrid,alen,1,ylen),1,1,1,2); % Terminal consumption
                        a1(:, T, :,:) = 0.0; % No assets after final period
                        v1(:, T, :,:) = model.utility(c1(:, T, :,:), par); % Terminal value function
                        alpha1(:, T,:,:) = 0.0;
                    else
                        for i = 1:ylen
                            yt = (T - age + 1 >= tr) * kappa * ygrid(i) + (T - age + 1 < tr) * ygrid(i);
                            for p = 1:alen
                                R=repmat((alpha*(1+r))+(1-alpha)*(1+r_values(rr)),alen,1);
                                ct = agrid(p) + yt - (amat ./ R); % Consumption equation
                                ct(agrid(p) < -par.borrow_limit) = par.borrow_limit; % Enforce borrowing constraint

                                ev = 0.5*(squeeze(v1(:,T-age+2,:,1))*pmat(i,:)') + 0.5*(squeeze(v1(:,T-age+2,:,2)) * pmat(i,:)');

                                vall = model.utility(ct, par) + beta * ev;
                                vall(ct <= 0.0) = -inf; % Penalize excessive borrowing
                                [vmax, ind] = max(vall,[],"all","linear");
                                [aind,alphaind] = ind2sub([alen 100],ind);

                               %  Policy functions (c1, a1) are computed by maximizing expected lifetime utility.
                                v1(p, T - age + 1, i,rr) = vmax;
                                c1(p, T - age + 1, i,rr) = ct(aind,alphaind); 
                                a1(p, T - age + 1, i,rr) = agrid(aind);
                                alpha1(p, T- age + 1, i, rr) = alpha(alphaind);
                            end
                        end
                    end
                end
            end
            sol.c = c1;
            sol.a = a1;
            sol.alpha = alpha1;
            sol.v = v1;
        end
    end
end