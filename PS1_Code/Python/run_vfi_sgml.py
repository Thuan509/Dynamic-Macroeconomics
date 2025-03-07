"""

run_vfi_dgm.py
--------------
This code solves the stochastic growth model using value function iteration.

"""

#%% Import from Python and set project directory
import os

# Set working directory
os.chdir(r"C:\Users\Do Thu An\OneDrive\Desktop\Dynamic Macroeconomics\Problem sets\Dynamic-Macroeconomics\PS1_Code\Python")
main = os.getcwd()
figout = os.path.join(main, "output", "figures")

#%% Import model functions
from model import planner
from solve import plan_allocations
from simulate import grow_economy
from my_graph import track_growth

#%% Stochastic Growth Model - Solving, Simulating, and Plotting
if __name__ == "__main__":
    print("\n-------------------- Running the Stochastic Growth Model --------------------\n")
    
    # Initialize model
    benevolent_dictator = planner()

    # Set parameters
    benevolent_dictator.setup(main=main, figout=figout)

    # Solve the model
    print("\nSolving the model...")
    plan_allocations(benevolent_dictator)

    # Simulate the model
    print("\nSimulating the economy...")
    grow_economy(benevolent_dictator)

    # Generate graphs
    print("\nGenerating graphs...")
    track_growth(benevolent_dictator)

    print("\n-------------------- Stochastic Growth Model Completed --------------------\n")