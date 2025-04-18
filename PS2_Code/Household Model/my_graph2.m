%% File Info.

%{

    my_graph.m
    ----------
    This code plots the value and policy functions and the time path of the variables.

%}

%% Graph class.

classdef my_graph2
    methods(Static)
        %% Plot value and policy functions.
        
        function [] = plot_policy(par,sol,sim, figout)
            %% Plot consumption policy function.

            ystate = par.ygrid;
            age = linspace(1,par.T,par.T);
            
            figure(1)
            
            surf(age(1:5:end),ystate,squeeze(sol.c(1,1:5:end,:))')
                xlabel({'t'},'Interpreter','latex')
                ylabel({'$y_{t}$'},'Interpreter','latex') 
                zlabel({'$c_{t}$'},'Interpreter','latex') 
            title('Consumption Policy Function, Lowest $a_t$','Interpreter','latex')
            saveas(gcf, fullfile(figout, 'cpol2_1.png'))

            
            figure(2)
            
            surf(age(1:5:end),ystate,squeeze(sol.c(end,1:5:end,:))')
                xlabel({'t'},'Interpreter','latex')
                ylabel({'$y_{t}$'},'Interpreter','latex') 
                zlabel({'$c_{t}$'},'Interpreter','latex') 
            title('Consumption Policy Function, Highest $a_t$','Interpreter','latex')
            saveas(gcf, fullfile(figout, 'cpol2_2.png'))
            
            %% Plot saving policy function.
            
            figure(3)
            
            surf(age(1:5:end),ystate,squeeze(sol.a(1,1:5:end,:))')
                xlabel({'t'},'Interpreter','latex')
                ylabel({'$y_{t}$'},'Interpreter','latex') 
                zlabel({'$a_{t+1}$'},'Interpreter','latex') 
            title('Saving Policy Function, Lowest $a_t$','Interpreter','latex')
            saveas(gcf, fullfile(figout, 'apol2_1.png'))

            
            figure(4)
            
            surf(age(1:5:end),ystate,squeeze(sol.a(end,1:5:end,:))')
                xlabel({'t'},'Interpreter','latex')
                ylabel({'$y_{t}$'},'Interpreter','latex') 
                zlabel({'$a_{t+1}$'},'Interpreter','latex') 
            title('Saving Policy Function, Highest $a_t$','Interpreter','latex')
            saveas(gcf, fullfile(figout, 'apol2_2.png'))
  

             %% Plot labor supply choice policy function.
            
            figure(5)
            
            surf(age(1:5:end),ystate,squeeze(sol.n(1,1:5:end,:))')
                xlabel({'t'},'Interpreter','latex')
                ylabel({'$y_{t}$'},'Interpreter','latex') 
                zlabel({'$n_{t+1}$'},'Interpreter','latex') 
            title('Labor Supply Choice Policy Function, Lowest $a_t$','Interpreter','latex')
            saveas(gcf, fullfile(figout, 'npol2_1.png'))
            
            figure(6)
            
            surf(age(1:5:end),ystate,squeeze(sol.n(end,1:5:end,:))')
                xlabel({'t'},'Interpreter','latex')
                ylabel({'$y_{t}$'},'Interpreter','latex') 
                zlabel({'$n_{t+1}$'},'Interpreter','latex') 
            title('Labor Supply Choice Policy Function, Highest $a_t$','Interpreter','latex')
            saveas(gcf, fullfile(figout, 'npol2_2.png'))
            
            %% Plot value function.
            
            figure(7)
            
            surf(age(1:5:end),ystate,squeeze(sol.v(1,1:5:end,:))')
                xlabel({'t'},'Interpreter','latex')
                ylabel({'$y_{t}$'},'Interpreter','latex') 
                zlabel({'$v_t(a_t,t)$'},'Interpreter','latex')
            title('Value Function, Lowest $a_t$','Interpreter','latex')
            saveas(gcf, fullfile(figout, 'vpol2_1.png'))

            figure(8)
            
            surf(age(1:5:end),ystate,squeeze(sol.v(end,1:5:end,:))')
                xlabel({'t'},'Interpreter','latex')
                ylabel({'$y_{t}$'},'Interpreter','latex') 
                zlabel({'$v_t(a_t,t)$'},'Interpreter','latex')
            title('Value Function, Highest $a_t$','Interpreter','latex')
            saveas(gcf, fullfile(figout, 'vpol2_2.png'))

            %% Plot consumption policy function.

            lcp_c = nan(par.T,1);
            lcp_a = nan(par.T,1);
            lcp_u = nan(par.T,1);
            lcp_n = nan(par.T,1);


            for i=1:par.T
                lcp_c(i) = mean(sim.csim(sim.tsim==i),"omitnan");
                lcp_a(i) = mean(sim.asim(sim.tsim==i),"omitnan");
                lcp_u(i) = mean(sim.usim(sim.tsim==i),"omitnan");
                lcp_n(i) = mean(sim.nsim(sim.tsim==i),"omitnan");
            end

            figure(9)
            
            plot(age,lcp_c)
                xlabel({'$Age$'},'Interpreter','latex')
                ylabel({'$c^{sim}_{t}$'},'Interpreter','latex') 
            title('LCP of Consumption')
            saveas(gcf, fullfile(figout, 'lcp_c2.png'))
            
            %% Plot saving policy function.
            
            figure(10)
            
            plot(age,lcp_a)
                xlabel({'$Age$'},'Interpreter','latex')
                ylabel({'$a^{sim}_{t+1}$'},'Interpreter','latex') 
            title('LCP of Savings')
            saveas(gcf, fullfile(figout, 'lcp_a2.png'))

            %% Plot labor supply choice policy function.
            
            figure(11)
            
            plot(age,lcp_n)
                xlabel({'$Age$'},'Interpreter','latex')
                ylabel({'$n^{sim}_{t+1}$'},'Interpreter','latex') 
            title('LCP of Labor Supply Choice')
            saveas(gcf, fullfile(figout, 'lcp_n2.png'))
            
            %% Plot value function.
            
            figure(12)
            
            plot(age,lcp_u)
                xlabel({'$Age$'},'Interpreter','latex')
                ylabel({'$u^{sim}_t$'},'Interpreter','latex') 
            title('LCP of Utility')
            saveas(gcf, fullfile(figout, 'lcp_u2.png'))

        end
        
    end
end
