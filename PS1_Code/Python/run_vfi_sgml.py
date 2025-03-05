"""

run_vfi_dgm.py
--------------
This code solves the stochastic growth model using value function iteration.

"""

#%% Import from Python and set project directory
import os
os.chdir("C:\\Users\\xmgb\\Dropbox\\02_FUV\\teaching\\spring_2025\\dynamic_macro\\code\\vfi_sgml_python")
main = os.getcwd()
figout = main+"\\output\\figures"

#%% Import from folder
from model import planner
from solve import plan_allocations
from simulate import grow_economy
from my_graph import track_growth

#%% Stochastic Growth Model.
benevolent_dictator = planner()

# Set the parameters, state space, and utility function.
benevolent_dictator.setup(main=main,figout=figout,beta = 0.96,sigma=2.00) # You can set the parameters here or use the defaults.

# Solve the model.
plan_allocations(benevolent_dictator) # Obtain the policy functions for capital.

# Simulate the model.
grow_economy(benevolent_dictator) # Simulate forward in time.

# Graphs.
track_growth(benevolent_dictator) # Plot policy functions and simulations.
