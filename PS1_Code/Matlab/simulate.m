%% File Info.
%{
    simulate.m
    ----------
    This code simulates the model.
%}

%% Simulate class.
classdef simulate
    methods(Static)
        %% Simulate the model.         
        function sim = grow(par, sol)            
            %% Set up.
            kgrid = par.kgrid; % Capital grid.
            Agrid = par.Agrid; % Productivity grid.
            yout = sol.y; % Production function.
            kpol = sol.k; % Capital policy function.
            cpol = sol.c; % Consumption policy function.
            ipol = sol.i; % Investment policy function.

            T = par.T; % Time periods.
            Asim = zeros(T*2,1); % Simulated productivity.
            ysim = zeros(T*2,1); % Simulated output.
            ksim = zeros(T*2,1); % Simulated capital.
            csim = zeros(T*2,1); % Simulated consumption.
            isim = zeros(T*2,1); % Simulated investment.
            usim = zeros(T*2,1); % Simulated utility.
            gsim = zeros(T*2,1); % Simulated government spending.
            tauksim = zeros(T*2,1); % Simulated capital tax.
            taunsim = zeros(T*2,1); % Simulated income tax.
            tsim = zeros(T*2,1); % Simulated transfers.

            %% Initialize random seed.
            rng(par.seed);

            % Stationary distribution for initial A.
            pmat0 = par.pmat^1000;
            pmat0 = pmat0(1,:);
            cmat = cumsum(par.pmat,2); % CDF matrix.

            % Draw initial conditions.
            k0_ind = randsample(par.klen,1); % Initial capital index.
            A0_ind = randsample(par.Alen,1,true,pmat0); % Initial productivity index.

            % Initial state.
            Asim(1) = Agrid(A0_ind);
            ksim(1) = kgrid(k0_ind);
            ksim2 = kpol(k0_ind, A0_ind);
            ysim(1) = yout(k0_ind,A0_ind);
            %isim(1) = ksim(1) * par.delta; % Investment follows capital accumulation
            isim(1) = ksim2 - (1 - par.delta) * ksim(1);
            gsim(1) = par.tauk * par.r * ksim(1) + par.taun * par.w * par.n - par.delta * par.tauk * ksim(1);
            tauksim(1) = par.tauk;
            taunsim(1) = par.taun;
            tsim(1) = par.t; % Government transfer
            csim(1) = cpol(k0_ind,A0_ind);
            usim(1) = model.utility(csim(1), gsim(1), par); % Utility function

            %% Simulate over time.
            for j = 2:T*2
                % Find closest index in kgrid.
                [~, kt_ind] = min(abs(kgrid - ksim(j-1))); % Ensure ksim updates correctly

                % Update states.
                Asim(j) = Agrid(A0_ind);
                %ksim(j) = kpol(kt_ind, A0_ind);
                ksim(j) = interp1(kgrid, kpol(:, A0_ind), ksim(j-1), 'linear');
                csim(j) = interp1(kgrid, cpol(:, A0_ind), ksim(j-1), 'linear');
                isim(j) = ksim(j) - (1 - par.delta) * ksim(j-1);
                ysim(j) = yout(kt_ind,A0_ind);
                %isim(j) = ksim(j) - (1 - par.delta) * ksim(j-1); % Corrected Investment Equation
                gsim(j) = par.tauk * par.r * ksim(j) + par.taun * par.w * par.n - par.delta * par.tauk * ksim(j);
                tauksim(j) = par.tauk;
                taunsim(j) = par.taun;
                tsim(j) = par.t;
                %csim(j) = cpol(kt_ind,A0_ind);
                usim(j) = model.utility(csim(j), gsim(j), par); % Fixed Utility

                % Draw next productivity state.
                A1_ind = find(rand <= cmat(A0_ind,:), 1, 'first');
                A0_ind = A1_ind;
            end

            %% Burn-in period.
            sim = struct();
            sim.Asim = Asim(T+1:2*T,1);
            sim.ysim = ysim(T+1:2*T,1);
            sim.ksim = ksim(T+1:2*T,1);
            sim.csim = csim(T+1:2*T,1);
            sim.isim = isim(T+1:2*T,1);
            sim.usim = usim(T+1:2*T,1);
            sim.gsim = gsim(T+1:2*T,1);
            sim.tauksim = tauksim(T+1:2*T,1); 
            sim.taunsim = taunsim(T+1:2*T,1); 
        end
    end
end
