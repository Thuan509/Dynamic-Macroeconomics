classdef simulate1
    methods(Static)
        function sim = firm_dynamics(par, sol)
            
            %% Set up

            kgrid = par.kgrid;
            Agrid = par.Agrid;
            pgrid = par.pgrid;

            vpol = sol.v;
            kpol = sol.k;
            rpol = sol.r;
            epol = sol.e;
            ppol = sol.p;
            xpol = sol.x;

            T = par.T;
            nlarge = par.nlarge;

            % Containers for simulations
            Asim = zeros(nlarge, T*2);
            ptsim = zeros(nlarge, T*2);
            vsim = zeros(nlarge, T*2);
            ksim = zeros(nlarge, T*2);
            esim = zeros(nlarge, T*2);
            psim = zeros(nlarge, T*2);
            xsim = zeros(nlarge, T*2);
            rsim = zeros(nlarge, T*2);

            %% Initial values
            rng(par.seed);
            pmat0 = par.pmat^1000;
            cmat = cumsum(par.pmat, 2);
            mtrix0 = par.mtrix^1000;
            cmtrix = cumsum(par.mtrix, 2);

            k0_ind = ones(nlarge, 1); 
            A0_ind = randsample(par.Alen, nlarge, true, pmat0(1,:));
            p0_ind = randsample(par.plen, nlarge, true, mtrix0(1,:));

            for i = 1:nlarge
                Asim(i,1) = Agrid(A0_ind(i));
                ptsim(i,1) = pgrid(p0_ind(i));
                vsim(i,1) = vpol(k0_ind(i), A0_ind(i), p0_ind(i));
                ksim(i,1) = kpol(k0_ind(i), A0_ind(i), p0_ind(i));
                xsim(i,1) = xpol(k0_ind(i), A0_ind(i), p0_ind(i));
                rsim(i,1) = rpol(k0_ind(i), A0_ind(i), p0_ind(i));
                esim(i,1) = epol(k0_ind(i), A0_ind(i), p0_ind(i));
                psim(i,1) = ppol(k0_ind(i), A0_ind(i), p0_ind(i));
            end

            %% Simulation
            for t = 2:T*2
                for i = 1:nlarge
                    kt_ind = find(kgrid == ksim(i, t-1));
                    kt_ind = kt_ind(1);
                    %if isempty(kt_ind), kt_ind = 1; end

                    Asim(i,t) = Agrid(A0_ind(i));
                    ptsim(i,t) = pgrid(p0_ind(i));
                    vsim(i,t) = vpol(kt_ind, A0_ind(i), p0_ind(i));
                    ksim(i,t) = kpol(kt_ind, A0_ind(i), p0_ind(i));
                    rsim(i,t) = rpol(kt_ind, A0_ind(i), p0_ind(i));
                    esim(i,t) = epol(kt_ind, A0_ind(i), p0_ind(i));
                    psim(i,t) = ppol(kt_ind, A0_ind(i), p0_ind(i));
                    xsim(i,t) = xpol(kt_ind, A0_ind(i), p0_ind(i));

                    % State transitions
                    A0_ind(i) = find(rand <= cmat(A0_ind(i), :), 1);
                    p0_ind(i) = find(rand <= cmtrix(p0_ind(i), :), 1);
                end
            end

            %% Store results (burn-in period)
            sim = struct();
            sim.Asim = Asim(:, T+1:end);
            sim.ptsim = ptsim(:, T+1:end);
            sim.vsim = vsim(:, T+1:end);
            sim.ksim = ksim(:, T+1:end);
            sim.rsim = rsim(:, T+1:end);
            sim.esim = esim(:, T+1:end);
            sim.psim = psim(:, T+1:end);
            sim.xsim = xsim(:, T+1:end);
        end
    end
end
