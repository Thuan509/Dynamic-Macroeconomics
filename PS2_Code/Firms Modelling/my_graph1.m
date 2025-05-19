%% File Info.
%{
    my_graph1.m
    ----------
    This code plots and saves the simulated averages and policy surfaces.
%}

%% Graph class.

classdef my_graph1
    methods(Static)

        %% 1. Plot and save simulated averages
        function [] = plot_policy(par, sim, figout, prefix)

            %% Create subfolder if it doesn't exist
            output_dir = fullfile(figout, prefix);
            if ~exist(output_dir, 'dir')
                mkdir(output_dir);
            end

            %% --- Averaged simulated time series ---

            % Time grid
            tgrid =  linspace(1,par.T,par.T);

            mean_Asim = nan(par.T, 1);
            mean_ksim = nan(par.T, 1);
            mean_rsim = nan(par.T, 1);
            mean_esim = nan(par.T, 1);
            mean_psim = nan(par.T, 1);
            mean_xsim = nan(par.T, 1);
            mean_ptsim = nan(par.T, 1);
            mean_vsim = nan(par.T, 1);

            for t = 1:par.T
                mean_Asim(t) = mean(sim.Asim(:, t), 'omitnan');
                mean_ksim(t) = mean(sim.ksim(:, t), 'omitnan');
                mean_rsim(t) = mean(sim.rsim(:, t), 'omitnan');
                mean_esim(t) = mean(sim.esim(:, t), 'omitnan');
                mean_psim(t) = mean(sim.psim(:, t), 'omitnan');
                mean_xsim(t) = mean(sim.xsim(:, t), 'omitnan');
                mean_ptsim(t) = mean(sim.ptsim(:, t), 'omitnan');
                mean_vsim(t) = mean(sim.vsim(:, t), 'omitnan');
            end

            %% --- Plot and Save ---

            figure(1)
            plot(tgrid, mean_Asim, 'LineWidth', 0.7)
            xlabel({'Time'},'Interpreter','latex')
            ylabel({'$\\bar{A}_t$'},'Interpreter','latex') 
            title('Average Simulated Revenue Shocks')
            saveas(gcf, fullfile(output_dir, 'simulated_mean_Asim.png'));

            figure(2)
            plot(tgrid, mean_ksim, 'LineWidth', 0.7)
            xlabel({'Time'},'Interpreter','latex')
            ylabel({'$\\bar{k}_t$'},'Interpreter','latex') 
            title('Average Simulated Capital Choice')
            saveas(gcf, fullfile(output_dir, 'simulated_mean_ksim.png'));

            figure(3)
            plot(tgrid, mean_rsim, 'LineWidth', 0.7)
            xlabel({'Time'},'Interpreter','latex')
            ylabel({'$\\bar{y}_t$'},'Interpreter','latex') 
            title('Average Simulated Revenue')
            saveas(gcf, fullfile(output_dir, 'simulated_mean_rsim.png'));

            figure(4)
            plot(tgrid, mean_psim, 'LineWidth', 0.7)
            xlabel({'Time'},'Interpreter','latex')
            ylabel({'$\\bar{\\pi}_t$'},'Interpreter','latex') 
            title('Average Simulated Profit')
            saveas(gcf, fullfile(output_dir, 'simulated_mean_psim.png'));

            figure(5)
            plot(tgrid, mean_vsim, 'LineWidth', 0.7)
            xlabel({'Time'},'Interpreter','latex')
            ylabel({'$\\bar{v}_t$'},'Interpreter','latex') 
            title('Average Simulated Firm Value')
            saveas(gcf, fullfile(output_dir, 'simulated_mean_vsim.png'));

            figure(6)
            plot(tgrid, mean_ptsim, 'LineWidth', 0.7)
            xlabel({'Time'},'Interpreter','latex')
            ylabel({'$\\bar{p}_t$'},'Interpreter','latex') 
            title('Average Simulated Price Choice')
            saveas(gcf, fullfile(output_dir, 'simulated_mean_ptsim.png'));

            figure(7)
            plot(tgrid, mean_xsim, 'LineWidth', 0.7)
            xlabel({'Time'},'Interpreter','latex')
            ylabel({'$\\bar{x}_t$'},'Interpreter','latex') 
            title('Average Simulated Input Cost')
            saveas(gcf, fullfile(output_dir, 'simulated_mean_xsim.png'));

            %% Plot simulated investment expenditure.

            figure(8)
            plot(tgrid, mean_esim, 'LineWidth', 0.7)
            xlabel({'Time'},'Interpreter','latex')
            ylabel({'$C(k^{sim}_{t+1},A^{sim}_t,k^{sim}_t)+\\pi^{sim}_t$'},'Interpreter','latex') 
            title('Simulated Investment Expenditure')
            saveas(gcf, fullfile(output_dir, 'simulated_investment_expenditure.png'));
        end

        %% 2. Plot and save 3D policy function surfaces
        function [] = plot_policy_3D(par, sol, figout, prefix)
            output_dir = fullfile(figout, prefix);
            if ~exist(output_dir, 'dir')
                mkdir(output_dir);
            end

            kgrid = par.kgrid;
            pgrid = par.pgrid;
            [K, P] = meshgrid(kgrid, pgrid);

            func_names = {'Value Function', 'Capital Policy', 'Investment Policy', ...
                          'Revenue', 'Expenditure', 'Profit', 'Input Cost'};
            func_fields = {'v', 'k', 'i', 'r', 'e', 'p', 'x'};

            for f = 1:length(func_fields)
                pol_func = sol.(func_fields{f});
                for a = [1, par.Alen]
                    surf_data = nan(length(pgrid), length(kgrid));
                    for i = 1:length(pgrid)
                        for j = 1:length(kgrid)
                            surf_data(i, j) = pol_func(j, a, i);
                        end
                    end

                    figure;
                    surf(K, P, surf_data);
                    xlabel('$k_t$', 'Interpreter', 'latex');
                    ylabel('$p_t$', 'Interpreter', 'latex');
                    zlabel(func_names{f}, 'Interpreter', 'latex');
                    title([func_names{f} ' at $A_t$ = ' num2str(par.Agrid(a))], 'Interpreter', 'latex');
                    grid on;
                    saveas(gcf, fullfile(output_dir, ['policy_' func_fields{f} '_A' num2str(a) '.png']));
                end
            end
        end

    end
end
