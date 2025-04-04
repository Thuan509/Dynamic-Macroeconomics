%% File Info.

%{
    simulate.m
    ----------
    This code simulates the model with multiple firms.
%}

%% Simulate class.

classdef simulate
    methods(Static)
        %% Simulate the model for multiple firms.
        
        function sim = firm_dynamics(par, sol)            
            %% Set up.
            
            kgrid = par.kgrid; % Capital today (state variable).
            Agrid = par.Agrid; % Productivity (state variable).

            vpol = sol.v; % Firm value.
            kpol = sol.k; % Policy function for capital.
            ipol = sol.i; % Policy function for investment.
            rpol = sol.r; % Optimal revenue.
            epol = sol.e; % Optimal total investment expenditure.
            ppol = sol.p; % Optimal profit.

            T = par.T; % Time periods.
            nfirm = par.nfirm; % Number of firms.

            % Containers for simulated values (each row is a firm, each column is a time step)
            Asim = zeros(nfirm, T);
            vsim = zeros(nfirm, T);
            ksim = zeros(nfirm, T);
            isim = zeros(nfirm, T);
            rsim = zeros(nfirm, T);
            esim = zeros(nfirm, T);
            psim = zeros(nfirm, T);
            
            %% Begin simulation.
            
            rng(par.seed);

            pmat0 = par.pmat^1000;
            pmat0 = pmat0(1,:); % Stationary distribution.
            cmat = cumsum(par.pmat, 2); % CDF matrix.

            % Initial conditions for all firms
            k0_ind = randsample(par.klen, nfirm, true); % Random initial capital stock indices
            A0_ind = randsample(par.Alen, nfirm, true, pmat0); % Random initial productivity indices

            % Initialize values at t=1
            for i = 1:nfirm
                Asim(i, 1) = Agrid(A0_ind(i)); % Initial productivity
                vsim(i, 1) = vpol(k0_ind(i), A0_ind(i)); % Initial firm value
                ksim(i, 1) = kpol(k0_ind(i), A0_ind(i)); % Initial capital stock
                isim(i, 1) = ipol(k0_ind(i), A0_ind(i)); % Initial investment
                rsim(i, 1) = rpol(k0_ind(i), A0_ind(i)); % Initial revenue
                esim(i, 1) = epol(k0_ind(i), A0_ind(i)); % Initial investment expenditure
                psim(i, 1) = ppol(k0_ind(i), A0_ind(i)); % Initial profit
            end

            %% Simulate endogenous and exogenous variables.

            for t = 2:T % Loop over time
                for i = 1:nfirm % Loop over firms
                    kt_ind = find(ksim(i, t-1) == kgrid); % Find capital index
                    Asim(i, t) = Agrid(A0_ind(i)); % Productivity in period t
                    vsim(i, t) = vpol(kt_ind, A0_ind(i)); % Firm value in period t
                    ksim(i, t) = kpol(kt_ind, A0_ind(i)); % Capital stock for period t+1
                    isim(i, t) = ipol(kt_ind, A0_ind(i)); % Investment in period t
                    rsim(i, t) = rpol(kt_ind, A0_ind(i)); % Revenue in period t
                    esim(i, t) = epol(kt_ind, A0_ind(i)); % Investment expenditure in period t
                    psim(i, t) = ppol(kt_ind, A0_ind(i)); % Profit in period t

                    % Draw next state for productivity
                    A1_ind = find(rand <= cmat(A0_ind(i), :));
                    A0_ind(i) = A1_ind(1); % Update state for next period
                end
            end

            %% Store results.
            sim = struct();
            sim.Asim = Asim; % Simulated productivity
            sim.vsim = vsim; % Simulated firm value
            sim.ksim = ksim; % Simulated capital choice
            sim.isim = isim; % Simulated investment
            sim.rsim = rsim; % Simulated revenue
            sim.esim = esim; % Simulated investment expenditure
            sim.psim = psim; % Simulated profit
        end
    end
end
