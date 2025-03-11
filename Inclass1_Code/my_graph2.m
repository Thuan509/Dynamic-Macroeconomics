%% File Info.

%{
    my_graph.m
    ----------
    This code plots policy functions and lifecycle simulation for the stochastic life cycle model.
%}

%% Graph class.

classdef my_graph2
    methods(Static)
        %% Plot policy functions by income state.
        
        function [] = plot_policy(par, sol)
            %% Plot consumption policy function by income state for different ages
            
            % Select sample periods to plot (working and retirement)
            periods_to_plot = [1, max(1, par.tr-1), par.tr, par.T];
            periods_names = {'Early Working', 'Late Working', 'Early Retirement', 'Late Retirement'};
            
            % Create figure with subplots
            figure('Position', [100, 100, 1200, 800])
            
            % Plot consumption policy
            for p = 1:length(periods_to_plot)
                t = periods_to_plot(p);
                
                % Create subplot
                subplot(2, 2, p)
                hold on
                
                % Plot for each income state
                for i_y = 1:par.ylen
                    plot(par.agrid, sol.c(:, i_y, t), 'LineWidth', 1.5)
                end
                
                hold off
                xlabel('Assets (a)')
                ylabel('Consumption (c)')
                title([periods_names{p}, ' (t=', num2str(t), ')'])
                
                if p == 1
                    legend_strs = cell(par.ylen, 1);
                    for i = 1:par.ylen
                        legend_strs{i} = ['y = ', num2str(par.ygrid(i), '%.2f')];
                    end
                    legend(legend_strs, 'Location', 'best')
                end
                
                grid on
            end
            
            sgtitle('Consumption Policy Function by Income State and Age')
            
            %% Plot savings policy function by income state for different ages
            
            figure('Position', [100, 100, 1200, 800])
            
            for p = 1:length(periods_to_plot)
                t = periods_to_plot(p);
                
                % Create subplot
                subplot(2, 2, p)
                hold on
                
                % Plot for each income state
                for i_y = 1:par.ylen
                    % Plot a' - a (savings)
                    plot(par.agrid, sol.a(:, i_y, t) - par.agrid, 'LineWidth', 1.5)
                end
                
                % Add zero line
                plot(par.agrid, zeros(size(par.agrid)), 'k--')
                
                hold off
                xlabel('Assets (a)')
                ylabel('Savings (a'' - a)')
                title([periods_names{p}, ' (t=', num2str(t), ')'])
                
                if p == 1
                    legend_strs = cell(par.ylen + 1, 1);
                    for i = 1:par.ylen
                        legend_strs{i} = ['y = ', num2str(par.ygrid(i), '%.2f')];
                    end
                    legend_strs{end} = 'Zero Savings';
                    legend(legend_strs, 'Location', 'best')
                end
                
                grid on
            end
            
            sgtitle('Savings Policy Function by Income State and Age')
        end
        
        %% Plot lifecycle simulation results.
        
        function [] = plot_lifecycle(par, sim)
            % Create time periods array (1 to T)
            periods = 1:par.T;
            retirement_line = par.tr;
            
            %% Plot all variables over the lifecycle
            figure('Position', [100, 100, 1000, 800])
            
            % Plot consumption
            subplot(4, 1, 1)
            plot(periods, sim.c, '-o', 'LineWidth', 1.5)
            xlabel('Age')
            ylabel('Consumption')
            title('Consumption over the Life Cycle')
            grid on
            % Add vertical line at retirement
            hold on
            xline(retirement_line, '--r', 'Retirement')
            hold off
            
            % Plot asset holdings
            subplot(4, 1, 2)
            plot(periods, sim.a, '-s', 'LineWidth', 1.5)
            xlabel('Age')
            ylabel('Assets')
            title('Asset Holdings over the Life Cycle')
            grid on
            % Add vertical line at retirement
            hold on
            xline(retirement_line, '--r', 'Retirement')
            hold off
            
            % Plot income
            subplot(4, 1, 3)
            plot(periods, sim.y, '-d', 'LineWidth', 1.5)
            xlabel('Age')
            ylabel('Income')
            title('Income over the Life Cycle')
            grid on
            % Add vertical line at retirement
            hold on
            xline(retirement_line, '--r', 'Retirement')
            hold off
            
            % Plot income state index
            subplot(4, 1, 4)
            plot(periods, sim.s, '-x', 'LineWidth', 1.5)
            xlabel('Age')
            ylabel('Income State')
            title('Income State over the Life Cycle')
            grid on
            % Add vertical line at retirement
            hold on
            xline(retirement_line, '--r', 'Retirement')
            hold off
            
            % Overall title
            sgtitle('Stochastic Life Cycle Model Simulation')
        end
        
        %% Plot multiple lifecycle simulations
        
        function [] = plot_multi_lifecycle(par, sims, n_sims)
            figure('Position', [100, 100, 1000, 800])
            
            % Plot consumption
            subplot(3, 1, 1)
            hold on
            for i = 1:n_sims
                plot(1:par.T, sims{i}.c, 'LineWidth', 0.5)
            end
            hold off
            xlabel('Age')
            ylabel('Consumption')
            title('Consumption Paths (Multiple Simulations)')
            grid on
            % Add vertical line at retirement
            hold on
            xline(par.tr, '--r', 'Retirement')
            hold off
            
            % Plot assets
            subplot(3, 1, 2)
            hold on
            for i = 1:n_sims
                plot(1:par.T, sims{i}.a, 'LineWidth', 0.5)
            end
            hold off
            xlabel('Age')
            ylabel('Assets')
            title('Asset Paths (Multiple Simulations)')
            grid on
            % Add vertical line at retirement
            hold on
            xline(par.tr, '--r', 'Retirement')
            hold off
            
            % Plot income
            subplot(3, 1, 3)
            hold on
            for i = 1:n_sims
                plot(1:par.T, sims{i}.y, 'LineWidth', 0.5)
            end
            hold off
            xlabel('Age')
            ylabel('Income')
            title('Income Paths (Multiple Simulations)')
            grid on
            % Add vertical line at retirement
            hold on
            xline(par.tr, '--r', 'Retirement')
            hold off
            
            sgtitle(['Stochastic Life Cycle Model: ', num2str(n_sims), ' Simulations'])
        end
    end
end