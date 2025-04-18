%% File Info.

%{
    my_graph1.m
    ----------
    This code plots the value and policy functions and the time path of the variables.
%}

%% Graph class.

classdef my_graph1
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
            sum_Asim = sum(sim.Asim, 1);
            sum_ksim = sum(sim.ksim, 1);
            sum_isim = sum(sim.isim, 1);
            sum_esim = sum(sim.esim, 1);
            sum_rsim = sum(sim.rsim, 1);
            sum_psim = sum(sim.psim, 1);
            sum_vsim = sum(sim.vsim, 1);

            %% Time grid.
            tgrid = 1:par.T;

            %% Plot averaged simulated revenue shocks.

            figure(7)
            plot(tgrid, sum_Asim, 'LineWidth', 1.5)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$\bar{A}_t$'},'Interpreter','latex') 
            title('Average Simulated Revenue Shocks')

            %% Plot averaged simulated capital choice.

            figure(8)
            plot(tgrid, sum_ksim, 'LineWidth', 1.5)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$\bar{k}_t$'},'Interpreter','latex') 
            title('Average Simulated Capital Choice')

            %% Plot averaged simulated investment expenditure.

            figure(9)
            plot(tgrid, sum_esim, 'LineWidth', 1.5)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$\bar{C}(k_{t+1},A_t,k_t) + \bar{\pi}_t$'},'Interpreter','latex') 
            title('Average Simulated Investment Expenditure')

            %% Plot averaged simulated investment.

            figure(10)
            plot(tgrid, sum_isim, 'LineWidth', 1.5)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$\bar{i}_t$'},'Interpreter','latex') 
            title('Average Simulated Investment')

            %% Plot averaged simulated revenue.

            figure(11)
            plot(tgrid, sum_rsim, 'LineWidth', 1.5)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$\bar{y}_t$'},'Interpreter','latex') 
            title('Average Simulated Revenue')

            %% Plot averaged simulated profit.

            figure(12)
            plot(tgrid, sum_psim, 'LineWidth', 1.5)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$\bar{\pi}_t$'},'Interpreter','latex') 
            title('Average Simulated Profit')

            %% Plot averaged simulated value function.

            figure(13)
            plot(tgrid, sum_vsim, 'LineWidth', 1.5)
                xlabel({'Time'},'Interpreter','latex')
                ylabel({'$\bar{v}_t$'},'Interpreter','latex') 
            title('Average Simulated Firm Value')

        end
    end
end
