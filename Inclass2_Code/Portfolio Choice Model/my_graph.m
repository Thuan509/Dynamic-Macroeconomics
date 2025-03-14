%% File Info.

%{

    my_graph.m
    ----------
    This code plots the value and policy functions and the time path of the variables.

%}

%% Graph class.

classdef my_graph
    methods(Static)
        function [] = plot_policy(par,sol,sim)
            astate = par.agrid;
            age = linspace(1,par.T,par.T);

            %% Plot Consumption Policy Function (Surf)
            figure(1)
            surf(age, astate, sol.c(:,:,1))
            xlabel({'$t$'},'Interpreter','latex')
            ylabel({'$a_{t}$'},'Interpreter','latex')
            zlabel({'$c_{t}$'},'Interpreter','latex')
            title('Consumption Policy Function','Interpreter','latex')

           

            %% Plot Saving Policy Function (Surf)
            figure(2)
            surf(age, astate, sol.a(:,:,1))
            xlabel({'$t$'},'Interpreter','latex')
            ylabel({'$a_{t}$'},'Interpreter','latex')
            zlabel({'$a_{t+1}$'},'Interpreter','latex')
            title('Saving Policy Function','Interpreter','latex')

            %% Plot Value Function (Surf)
            figure(3)
            surf(age, astate, sol.v(:,:,1))
            xlabel({'$t$'},'Interpreter','latex')
            ylabel({'$a_{t}$'},'Interpreter','latex')
            zlabel({'$v_{t}$'},'Interpreter','latex')
            title('Value Function','Interpreter','latex')

            %% Plot Asset Portfolio Policy Function (Surf)
            figure(5)
            surf(age, astate, sol.alpha(:,:,1))
            xlabel({'$t$'},'Interpreter','latex')
            ylabel({'$a_{t}$'},'Interpreter','latex')
            zlabel({'$\alpha$'},'Interpreter','latex')
            title('Asset Portfolio Policy Function','Interpreter','latex')

            %% Life Cycle Profiles (Average across Simulations)
            lcp_c = mean(sim.c,2);
            lcp_a = mean(sim.a,2);
            lcp_y = mean(sim.y,2);
            lcp_r = mean(sim.r,2);

            figure(6)
            plot(age, lcp_c)
            xlabel({'Age'},'Interpreter','latex')
            ylabel({'$c_t$'},'Interpreter','latex')
            title('Life Cycle Profile of Consumption')

            figure(7)
            plot(age, lcp_a)
            xlabel({'Age'},'Interpreter','latex')
            ylabel({'$a_t$'},'Interpreter','latex')
            title('Life Cycle Profile of Assets')

            figure(8)
            plot(age, lcp_y)
            xlabel({'Age'},'Interpreter','latex')
            ylabel({'$y_t$'},'Interpreter','latex')
            title('Life Cycle Profile of Income')

            figure(9)
            plot(age, lcp_r)
            xlabel({'Time'}, 'Interpreter', 'latex')
            ylabel({'$r_t$'}, 'Interpreter', 'latex')
            title('Life Cycle Interest Rate')
        end
    end
end
