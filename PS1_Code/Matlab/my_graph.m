%% File Info.
%{
    my_graph.m
    ----------
    This code plots the value and policy functions.
%}

%% Graph class.

classdef my_graph
    methods(Static)
        %% Plot value and policy functions.
        
        function [] = plot_policy(par, sol, sim, figout)
            %% Plot production function.
            
            figure(1)
            plot(par.kgrid, sol.y)
            xlabel({'$k_{t}$'}, 'Interpreter', 'latex')
            ylabel({'$y_{t}$'}, 'Interpreter', 'latex')
            title('Production Function')
            savefig(strcat(figout, 'ypol.fig'))
            
            %% Plot capital policy function.
            
            figure(2)
            plot(par.kgrid, sol.k)
            xlabel({'$k_{t}$'}, 'Interpreter', 'latex')
            ylabel({'$k_{t+1}$'}, 'Interpreter', 'latex')
            title('Capital Policy Function')
            savefig(strcat(figout, 'kpol.fig'))
            
            %% Plot consumption policy function.
            
            figure(3)
            plot(par.kgrid, sol.c)
            xlabel({'$k_{t}$'}, 'Interpreter', 'latex')
            ylabel({'$c_{t}$'}, 'Interpreter', 'latex')
            title('Consumption Policy Function')
            savefig(strcat(figout, 'cpol.fig'))
            
            %% Plot investment policy function.
            
            figure(4)
            plot(par.kgrid, sol.i)
            xlabel({'$k_{t}$'}, 'Interpreter', 'latex')
            ylabel({'$i_{t}$'}, 'Interpreter', 'latex')
            title('Investment Policy Function')
            savefig(strcat(figout, 'ipol.fig'))
            
            %% Plot value function.
            
            figure(5)
            plot(par.kgrid, sol.v)
            xlabel({'$k_{t}$'}, 'Interpreter', 'latex')
            ylabel({'$v_t(k_t,A_t)$'}, 'Interpreter', 'latex')
            title('Value Function')
            savefig(strcat(figout, 'vfun.fig'))
            
            %% Define time grid.
            tgrid = 1:par.T;
            
            %% Plot simulated output.
            figure(6)
            plot(tgrid, sim.ysim)
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$y^{sim}_t$'}, 'Interpreter', 'latex')
            title('Simulated Output')
            savefig(strcat(figout, 'ysim.fig'))
            
            %% Plot simulated capital choice.
            figure(7)
            plot(tgrid, sim.ksim)
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$k^{sim}_t$'}, 'Interpreter', 'latex')
            title('Simulated Capital Choice')
            savefig(strcat(figout, 'ksim.fig'))
            
            %% Plot simulated consumption.
            figure(8)
            plot(tgrid, sim.csim)
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$c^{sim}_t$'}, 'Interpreter', 'latex')
            title('Simulated Consumption')
            savefig(strcat(figout, 'csim.fig'))
            
            %% Plot simulated investment.
            figure(9)
            plot(tgrid, sim.isim)
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$i^{sim}_t$'}, 'Interpreter', 'latex')
            title('Simulated Investment')
            savefig(strcat(figout, 'isim.fig'))
            
            %% Plot simulated utility.
            figure(10)
            plot(tgrid, sim.usim)
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$u^{sim}_t$'}, 'Interpreter', 'latex')
            title('Simulated Utility')
            savefig(strcat(figout, 'usim.fig'))
            
            %% Plot simulated productivity.
            figure(11)
            plot(tgrid, sim.Asim)
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$A^{sim}_t$'}, 'Interpreter', 'latex')
            title('Simulated Productivity')
            savefig(strcat(figout, 'Asim.fig'))
            
            %% Plot simulated government spending.
            figure(12)
            plot(tgrid, sim.gsim)
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$G^{sim}_t$'}, 'Interpreter', 'latex')
            title('Simulated Government Spending')
            savefig(strcat(figout, 'Gsim.fig'))
            
        end
    end
end