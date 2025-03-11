%% File Info.

%{
    simulate.m
    ----------
    This code simulates the life cycle model with stochastic income.
%}

%% Simulate class.

classdef simulate2
    methods(Static)
        %% Simulate the model.
        
        function sim = grow(par, sol)
            %% Set up.
            
            agrid = par.agrid;    % Asset grid
            ygrid = par.ygrid;    % Income grid
            ytrans = par.ytrans;  % Income transition matrix
            T = par.T;            % Total lifespan
            tr = par.tr;          % Retirement age
            r = par.r;            % Interest rate
            kappa = par.kappa;    % Pension fraction
            
            a_pol = sol.a;        % Policy function for assets
            c_pol = sol.c;        % Policy function for consumption
            
            % Containers for simulation
            a = zeros(T, 1);      % Assets by age
            c = zeros(T, 1);      % Consumption by age
            y = zeros(T, 1);      % Income by age
            s = zeros(T, 1);      % Income state index by age
            
            %% Initialize random number generator
            
            rng(par.seed);
            
            %% Initial conditions
            
            a(1) = par.a0;                            % Initial assets
            s(1) = randi(par.ylen);                   % Random initial income state
            y(1) = ygrid(s(1));                       % Initial income level
            
            %% Simulate lifecycle
            
            for t = 1:T
                % Interpolate to find optimal policy
                [a_idx, ~] = model2.interp_index(a(t), agrid);
                
                if t < tr
                    % Working period: stochastic income
                    [c(t), a_next] = model2.interp_policy(a(t), s(t), t, agrid, a_pol, c_pol);
                    
                    % Draw next period's income state if not moving to retirement
                    if t < T && t < tr-1
                        cdf = cumsum(ytrans(s(t),:));
                        s(t+1) = find(rand <= cdf, 1, 'first');
                        y(t+1) = ygrid(s(t+1));
                    elseif t < T && t == tr-1
                        % Last working period - pension will be based on this income
                        s(t+1) = s(t);  % Keep the state index (for consistency)
                        y(t+1) = kappa * y(t);  % Pension starts
                    end
                else
                    % Retirement period: fixed pension (already set in previous periods)
                    [c(t), a_next] = model2.interp_policy(a(t), s(t), t, agrid, a_pol, c_pol);
                    
                    if t < T
                        s(t+1) = s(t);  % Keep the state index (for consistency)
                        y(t+1) = y(t);  % Pension is constant during retirement
                    end
                end
                
                % Update assets for next period if not the last period
                if t < T
                    a(t+1) = a_next;
                end
            end
            
            %% Calculate utilities
            
            u = zeros(T, 1);
            for t = 1:T
                u(t) = model.utility(c(t), par);
            end
            
            %% Store results
            
            sim = struct();
            sim.a = a;      % Asset holdings
            sim.c = c;      % Consumption
            sim.y = y;      % Income
            sim.s = s;      % Income state index
            sim.u = u;      % Utility
        end
    end
end