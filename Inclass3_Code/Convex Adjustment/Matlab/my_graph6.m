%% File Info.

%{
    my_graph.m
    ----------
    This code plots the value and policy functions and the time path of the variables.
%}

%% Graph class.

classdef my_graph
    methods(Static)
        %% Plot value and policy functions.
        
        function [] = plot_policy(par,sol,sim)
            %% Plot capital policy function.

            figure(1)
            
            plot(par.kgrid,sol.k)
                xlabel({'$k_t$'},'Interpreter','latex')
                ylabel({'$k_{t+1}$'},'Interpreter','latex') 
            title('Capital Policy Function','Interpreter','latex')
            
            %% Plot investment policy function.
            
            figure(2)
            
            plot(par.kgrid,sol.i)
                xlabel({'$k_t$'},'Interpreter','latex')
                ylabel({'$i_{t}$'},'Interpreter','latex') 
            title('Investment Policy Function','Interpreter','latex')
            
            %% Plot revenue function.
            
            figure(3)
            
            plot(par.kgrid,sol.r)
                xlabel({'$k_t$'},'Interpreter','latex')
                ylabel({'$r_{t}$'},'Interpreter','latex') 
            title('Revenue Function','Interpreter','latex')
            
            %% Plot expenditure function.
            
            figure(4)
            
            plot(par.kgrid,sol.e)
                xlabel({'$k_t$'},'Interpreter','latex')
                ylabel({'$C(k_{t+1},A_t,k_t)+pI_t$'},'Interpreter','latex') 
            title('Expenditure Function','Interpreter','latex')
            
            %% Plot profit function.
            
            figure(5)
            
            plot(par.kgrid,sol.p)
                xlabel({'$k_t$'},'Interpreter','latex')
                ylabel({'$C(k_{t+1},A_t,k_t)+\pi_t$'},'Interpreter','latex') 
            title('Profit Function','Interpreter','latex')
            
            %% Plot value function.
            
            figure(6)
            
            plot(par.kgrid,sol.v)
                xlabel({'$k_t$'},'Interpreter','latex')
                ylabel({'$v_{t}$'},'Interpreter','latex') 
            title('Value Function','Interpreter','latex')
            
            %% Compute averages over firms
            avg_Asim = mean(sim.Asim, 1);
            avg_ksim = mean(sim.ksim, 1);
            avg_isim = mean(sim.isim, 1);
            avg_esim = mean(sim.esim, 1);
            avg_rsim = mean(sim.rsim, 1);
            avg_psim = mean(sim.psim, 1);
            avg_vsim = mean(sim.vsim, 1);

            %% Time grid.
            tgrid = 1:par.T;

            %% Plot averaged simulated revenue shocks.

            figure(7)
            plot(tgrid, avg_Asim, 'LineWidth', 1.5)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$\bar{A}_t$'},'Interpreter','latex') 
            title('Average Simulated Revenue Shocks')

            %% Plot averaged simulated capital choice.

            figure(8)
            plot(tgrid, avg_ksim, 'LineWidth', 1.5)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$\bar{k}_t$'},'Interpreter','latex') 
            title('Average Simulated Capital Choice')

            %% Plot averaged simulated investment expenditure.

            figure(9)
            plot(tgrid, avg_esim, 'LineWidth', 1.5)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$\bar{C}(k_{t+1},A_t,k_t) + \bar{\pi}_t$'},'Interpreter','latex') 
            title('Average Simulated Investment Expenditure')

            %% Plot averaged simulated investment.

            figure(10)
            plot(tgrid, avg_isim, 'LineWidth', 1.5)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$\bar{i}_t$'},'Interpreter','latex') 
            title('Average Simulated Investment')

            %% Plot averaged simulated revenue.

            figure(11)
            plot(tgrid, avg_rsim, 'LineWidth', 1.5)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$\bar{y}_t$'},'Interpreter','latex') 
            title('Average Simulated Revenue')

            %% Plot averaged simulated profit.

            figure(12)
            plot(tgrid, avg_psim, 'LineWidth', 1.5)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$\bar{\pi}_t$'},'Interpreter','latex') 
            title('Average Simulated Profit')

            %% Plot averaged simulated value function.

            figure(13)
            plot(tgrid, avg_vsim, 'LineWidth', 1.5)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$\bar{v}_t$'},'Interpreter','latex') 
            title('Average Simulated Firm Value')

        end
    end
end
