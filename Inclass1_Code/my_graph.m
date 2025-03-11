%% File Info.

%{
    my_graph.m
    ----------
    This code plots policy functions and lifecycle simulation for the life cycle model.
%}

%% Graph class.

classdef my_graph
    methods(Static)
        %% Plot policy functions.
        
        function [] = plot_policy(par, sol)
            %% Plot consumption function for different periods
            
            figure(1)
            hold on
            
            % Plot consumption policy for each period
            for t = 1:par.T
                if t < par.tr
                    plot(par.agrid, sol.c(:,t), 'LineWidth', 1.5)
                else
                    plot(par.agrid, sol.c(:,t), '--', 'LineWidth', 1.5)
                end
            end
            
            hold off
            xlabel('Assets (a_t)')
            ylabel('Consumption (c_t)') 
            title('Consumption Policy Function by Age')
            
            % Create legend based on working/retirement status
            legend_strs = cell(par.T, 1);
            for t = 1:par.T
                if t < par.tr
                    legend_strs{t} = ['Working Age ' num2str(t)];
                else
                    legend_strs{t} = ['Retirement Age ' num2str(t)];
                end
            end
            legend(legend_strs, 'Location', 'best')
            grid on
            
            %% Plot asset policy function for different periods
            
            figure(2)
            hold on
            
            % Plot asset policy for each period (except last one)
            for t = 1:(par.T-1)
                if t < par.tr
                    plot(par.agrid, sol.a(:,t), 'LineWidth', 1.5)
                else
                    plot(par.agrid, sol.a(:,t), '--', 'LineWidth', 1.5)
                end
            end
            
            % Add 45-degree line for reference
            plot(par.agrid, par.agrid, 'k:', 'LineWidth', 1.0)
            
            hold off
            xlabel('Current Assets (a_t)')
            ylabel('Next Period Assets (a_{t+1})') 
            title('Asset Policy Function by Age')
            
            % Create legend based on working/retirement status
            legend_strs = cell(par.T, 1);
            for t = 1:(par.T-1)
                if t < par.tr
                    legend_strs{t} = ['Working Age ' num2str(t)];
                else
                    legend_strs{t} = ['Retirement Age ' num2str(t)];
                end
            end
            legend_strs{par.T} = '45° Line';
            legend(legend_strs, 'Location', 'best')
            grid on
        end
        
        %% Plot lifecycle simulation results.
        
        function [] = plot_lifecycle(par, sim)
            % Create time periods array (1 to T)
            periods = 1:par.T;
            
            %% Plot all variables over the lifecycle
            figure(3)
            
            % Plot consumption
            subplot(3,1,1)
            plot(periods, sim.c, '-o', 'LineWidth', 2)
            xlabel('Age')
            ylabel('Consumption')
            title('Consumption over the Life Cycle')
            grid on
            % Add vertical line at retirement
            hold on
            xline(par.tr, '--r', 'Retirement')
            hold off
            
            % Plot asset holdings
            subplot(3,1,2)
            plot(periods, sim.a, '-s', 'LineWidth', 2)
            xlabel('Age')
            ylabel('Assets')
            title('Asset Holdings over the Life Cycle')
            grid on
            % Add vertical line at retirement
            hold on
            xline(par.tr, '--r', 'Retirement')
            hold off
            
            % Plot income
            subplot(3,1,3)
            plot(periods, sim.y, '-d', 'LineWidth', 2)
            xlabel('Age')
            ylabel('Income')
            title('Income over the Life Cycle')
            grid on
            % Add vertical line at retirement
            hold on
            xline(par.tr, '--r', 'Retirement')
            hold off
            
            % Overall title
            sgtitle('Life Cycle Model Simulation')
        end
    end
end