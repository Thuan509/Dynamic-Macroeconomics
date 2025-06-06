%% File Info.

%{

    simulate.m
    ----------
    This code simulates the model.

%}

%% Simulate class.

classdef simulate2
    methods(Static)
        %% Simulate the model. 
        
        function sim = lc(par,sol)            
            %% Set up.
            
            agrid = par.agrid; % Assets today (state variable).

            apol = sol.a; % Policy function for capital.
            cpol = sol.c; % Policy function for consumption.
            npol = sol.n; % Policy function for labor.
            n_all = sol.n_all;

            TT = par.TT; % Time periods.
            NN = par.NN; % People.
            T = par.T; % Life span.
            tr = par.tr; % Retirement.

            kappa = par.kappa; % Share of income as pension.
            ygrid = par.ygrid; % Exogenous income.
            pmat = par.pmat; % Transition matrix.
            Gmat = par.Gmat; % Matrix for G_t
            

            ysim = nan(TT,NN); % Container for simulated income.
            asim = nan(TT,NN); % Container for simulated savings.
            tsim = nan(TT,NN); % Container for simulated age.
            csim = nan(TT,NN); % Container for simulated consumption.
            usim = nan(TT,NN); % Container for simulated utility.
            nsim = nan(TT,NN); % Container for simulated labor choice.

            %% Begin simulation.
            
            rng(par.seed);

            pmat0 = pmat^100; % Stationary distribution.
            cmat = cumsum(pmat,2); % CDF matrix.

            y0_ind = randsample(par.ylen,NN,true,pmat0(1,:))'; % Index for initial income.
            a0_ind = randsample(par.alen,NN,true)'; % Index for initial wealth.
            t0_ind = ones(T,NN)'; % Index for initial wealth.
            yr = nan(NN,1); % Retirement income.

            for i = 1:NN % Person loop.
                n = npol(a0_ind(i),t0_ind(i),y0_ind(i));
                nsim(1,i) = npol(a0_ind(i),t0_ind(i),y0_ind(i));

                
                if t0_ind(i)>=tr % Retired now.
                    nt_T = n_all(y0_ind(i), a0_ind(i));
                    yr(i) = ygrid(y0_ind(i)); % Store for pension.
                    ysim(1,i) = kappa.*yr(i)*Gmat(tr-1)*nt_T; % Pension in period 0 given age.
                    %n = 0.0;
                else
                    ysim(1,i) = Gmat(1) * ygrid(y0_ind(i)) * n; % Pension in period 0 given age.
                end

                tsim(1,i) = t0_ind(i); % Age in period 0.
                csim(1,i) = cpol(a0_ind(i),t0_ind(i),y0_ind(i)); % Consumption in period 0 given a0.
                asim(1,i) = apol(a0_ind(i),t0_ind(i),y0_ind(i)); % Savings for period 1 given a0.
                %usim(1,i) = model.utility(csim(1,i),n,par); % Utility in period 0 given a0.
                %nsim(1,i) = npol(a0_ind(i),t0_ind(i),y0_ind(i)); % Savings for period 1 given a0.

                if t0_ind(i) == tr-1 % Retired next period.
                    yr(i) = ygrid(y0_ind(i)); % Store as pension for next period
                elseif t0_ind(i) < tr-1
                    y1_ind = find(rand<=cmat(y0_ind(i),:)); % Draw income shock for next period.
                    y0_ind(i) = y1_ind(1);
                end

            end

            usim(1, :) = model2.utility(csim(1, :), nsim(1,:), par);

            %% Simulate endogenous variables.

            for j = 2:TT % Time loop.
                for i = 1:NN % Person loop.

                    age = tsim(j-1,i)+1; % Age in period t.

                    if age <= T % Check if still alive.
                        tsim(j,i) = age; % Age in period t.
                        at_ind = find(asim(j-1,i)==agrid); % Savings choice in the previous period is the state today. Find where the latter is on the grid.
                        n = npol(at_ind,age,y0_ind(i));
                        nsim(j,i) = npol(at_ind,age,y0_ind(i));

                        if age>=tr % Retired
                            nt_T = n_all(y0_ind(i), a0_ind(i));
                            ysim(j,i) = kappa.*yr(i)*Gmat(tr-1)*nt_T; % Pension in period t given age.
                            %n=0.0;
                        else
                            ysim(j,i) = Gmat(age)*ygrid(y0_ind(i))*n; % Pension in period t given age.
                            
                        end

                        %[~, at_ind] = min(abs(agrid - asim(j-1,i)));
                        %at_ind = max(1, min(par.alen, at_ind)); % clamp index within [1, alen]
                        csim(j,i) = cpol(at_ind,age,y0_ind(i)); % Consumption in period t.
                        asim(j,i) = apol(at_ind,age,y0_ind(i)); % Savings for period t+1.
                        %nsim(j,i) = npol(at_ind,age,y0_ind(i));
                        usim(j,i) = model2.utility(csim(j,i),nsim(j,i),par); % Utility in period t.

                        if age == tr-1 % Retire next period
                            yr(i) = ygrid(y0_ind(i)); % Store as pension for next period
                        elseif age < tr-1
                            y1_ind = find(rand<=cmat(y0_ind(i),:),1); % Draw income shock for next period.
                            y0_ind(i) = y1_ind(1);
                        end

                    end
                end
            end

            sim = struct();
            
            sim.ysim = ysim; % Simulated output.
            sim.asim = asim; % Simulated savings.
            sim.tsim = tsim; % Simulated age.
            sim.csim = csim; % Simulated consumption.
            sim.usim = usim; % Simulated utility.
            sim.nsim = nsim;
             
        end
        
    end
end