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
            saveas(gcf, fullfile(figout, 'ypol.png'))
            
            %% Plot capital policy function.
            
            figure(2)
            plot(par.kgrid, sol.k)
            xlabel({'$k_{t}$'}, 'Interpreter', 'latex')
            ylabel({'$k_{t+1}$'}, 'Interpreter', 'latex')
            title('Capital Policy Function')
            saveas(gcf, fullfile(figout, 'kpol.png'))
            
            %% Plot consumption policy function.
            
            figure(3)
            plot(par.kgrid, sol.c)
            xlabel({'$k_{t}$'}, 'Interpreter', 'latex')
            ylabel({'$c_{t}$'}, 'Interpreter', 'latex')
            title('Consumption Policy Function')
            saveas(gcf, fullfile(figout, 'cpol.png'))
            
            %% Plot investment policy function.
            
            figure(4)
            plot(par.kgrid, sol.i)
            xlabel({'$k_{t}$'}, 'Interpreter', 'latex')
            ylabel({'$i_{t}$'}, 'Interpreter', 'latex')
            title('Investment Policy Function')
            saveas(gcf, fullfile(figout, 'ipol.png'))
            
            %% Plot value function.
            
            figure(5)
            plot(par.kgrid, sol.v)
            xlabel({'$k_{t}$'}, 'Interpreter', 'latex')
            ylabel({'$v_t(k_t,A_t)$'}, 'Interpreter', 'latex')
            title('Value Function')
            saveas(gcf, fullfile(figout, 'vfun.png'))

            %% Plot government policy function.
            
            figure(6)
            plot(par.kgrid, sol.g)
            xlabel({'$g_{t}$'}, 'Interpreter', 'latex')
            ylabel({'$g_{t+1}$'}, 'Interpreter', 'latex')
            title('Government Policy Function')
            saveas(gcf, fullfile(figout, 'gpol.png'))
            
            %% Define time grid.
            tgrid = 1:par.T;
            
            %% Plot simulated output.
            figure(7)
            plot(tgrid, sim.ysim)
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$y^{sim}_t$'}, 'Interpreter', 'latex')
            title('Simulated Output')
            saveas(gcf, fullfile(figout, 'ysim.png'))
            
            %% Plot simulated capital choice.
            figure(8)
            plot(tgrid, sim.ksim)
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$k^{sim}_t$'}, 'Interpreter', 'latex')
            title('Simulated Capital Choice')
            saveas(gcf, fullfile(figout, 'ksim.png'))
            
            %% Plot simulated consumption.
            figure(9)
            plot(tgrid, sim.csim)
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$c^{sim}_t$'}, 'Interpreter', 'latex')
            title('Simulated Consumption')
            saveas(gcf, fullfile(figout, 'csim.png'))
            
            %% Plot simulated investment.
            figure(10)
            plot(tgrid, sim.isim)
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$i^{sim}_t$'}, 'Interpreter', 'latex')
            title('Simulated Investment')
            saveas(gcf, fullfile(figout, 'isim.png'))
            
            %% Plot simulated utility.
            figure(11)
            plot(tgrid, sim.usim)
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$u^{sim}_t$'}, 'Interpreter', 'latex')
            title('Simulated Utility')
            saveas(gcf, fullfile(figout, 'usim.png'))
            
            %% Plot simulated productivity.
            figure(12)
            plot(tgrid, sim.Asim)
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$A^{sim}_t$'}, 'Interpreter', 'latex')
            title('Simulated Productivity')
            saveas(gcf, fullfile(figout, 'Asim.png'))
            
            %% Plot simulated government spending.
            figure(13)
            plot(tgrid, sim.gsim)
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$G^{sim}_t$'}, 'Interpreter', 'latex')
            title('Simulated Government Spending')
            saveas(gcf, fullfile(figout, 'Gsim.png'))
            
        end
    end
end