classdef realdata_plot
    methods(Static)
        function plot_economic_data()
            %% Extract Variables
            year = data.Year;
            gdp = data.GDP_constant_2015_USD_;
            gdp_growth = data.GDP_growth_rate__;
            capital_share = data.Gross_fixed_capital_formation__of_GDP_;
            capital_growth = data.Gross_capital_formation_growth__;
            fdi = data.FDI_inflow_BoP_US_;
            export_growth = data.Exports_of_goods_and_Services_annual__growth_;
            
            %% Create Output Folder
            output_folder = "Output";
            if ~exist(output_folder, 'dir')
                mkdir(output_folder);
            end
            
            %% Plot GDP
            figure;
            plot(year, gdp, 'b-', 'LineWidth', 1.5);
            title('GDP (Billion US$)');
            xlabel('Year'); 
            ylabel('Billion US ($)');
            grid on;
            saveas(gcf, fullfile(output_folder, 'gdp.png'));
            
            %% Plot GDP Growth Rate
            figure;
            plot(year, gdp_growth, 'b-', 'LineWidth', 1.5);
            title("Vietnam's GDP Growth Rate");
            xlabel('Year'); 
            ylabel('GDP Growth Rate (%)');
            grid on;
            saveas(gcf, fullfile(output_folder, 'gdp_growth_rate.png'));
            
            %% Plot Capital Share
            figure;
            plot(year, capital_share, 'r-', 'LineWidth', 1.5);
            yline(mean(capital_share), '--k'); % Average line
            title('Gross Fixed Capital Formation (% of GDP)');
            xlabel('Year'); 
            ylabel('% Share of GDP');
            grid on;
            saveas(gcf, fullfile(output_folder, 'capital_share.png'));
            
            %% Plot Capital Formation & GDP Growth
            figure;
            [ax, h1, h2] = plotyy(year, capital_growth, year, gdp_growth);
            title('Vietnam: Capital Formation & GDP Growth');
            xlabel('Year');
            ylabel(ax(1), 'Gross Capital Formation Growth (%)');
            ylabel(ax(2), 'GDP Growth Rate (%)');
            grid on;
            saveas(gcf, fullfile(output_folder, 'capital_growth.png'));
            
            %% Plot Incremental Capital-Output Ratio (ICOR)
            investment = (capital_share / 100) .* gdp;
            change_in_gdp = gdp(2:end) - gdp(1:end-1);
            icor = investment(2:end) ./ change_in_gdp;
            figure;
            plot(year(2:end), icor, 'g-', 'LineWidth', 1.5);
            title('Incremental Capital-Output Ratio (ICOR)');
            xlabel('Year'); 
            ylabel('ICOR');
            grid on;
            saveas(gcf, fullfile(output_folder, 'icor.png'));
            
            %% Plot Export Growth
            figure;
            plot(year, export_growth, 'm-', 'LineWidth', 1.5);
            title("Vietnam's Export Growth");
            xlabel('Year'); 
            ylabel('Exports (Annual % Growth)');
            grid on;
            saveas(gcf, fullfile(output_folder, 'export_growth.png'));
            
            %% Plot FDI Inflow
            figure;
            plot(year, fdi, 'g-', 'LineWidth', 1.5);
            title("Vietnam's FDI Inflow (BoP, Billion US$)");
            xlabel('Year'); 
            ylabel('FDI (Billion US$)');
            grid on;
            saveas(gcf, fullfile(output_folder, 'fdi.png'));
            
            fprintf("Plots saved in folder: %s\n", output_folder);
        end
    end
end