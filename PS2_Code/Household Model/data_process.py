# Import libraries
import pandas as pd
import numpy as np
import os

# Define main project folder
main = r'C:\Users\Do Thu An\OneDrive\Desktop\Dynamic Macroeconomics\Problem sets\Dynamic-Macroeconomics\PS2_Code'
# Set project folder as current working directory
os.chdir(main)

# Define data file path
data_path = os.path.join(main, 'Data Files', 'VHLSS 2008 Data')

# Load muc123a
muc123a = pd.read_csv(os.path.join(data_path, 'muc123a.csv'))

# Create household size column
muc123a['hsize'] = muc123a.groupby(['tinh', 'huyen', 'xa', 'diaban', 'hoso'])['matv'].transform('max')

# Filter to keep only rows where m1ac3 == 1 (household head)
# Keep only household heads who are male and age ≥ 25
muc123a = muc123a[(muc123a['m1ac3'] == 1) & (muc123a['m1ac2'] == 1) & (muc123a['m1ac5'] >= 25)]

# Only keep the variables of label and individuals age 
columns123a = ['tinh', 'huyen', 'xa', 'diaban', 'hoso', 'matv', 'hsize', 'm1ac2', 'm1ac3', 'm1ac5']

# We take m123a as the base to merge the labels of other files
# Select the columns from muc123a to create a new dataframe
df = muc123a[columns123a].copy()

# ************************ CALCULATE TOTAL HOUSEHOLD INCOME *******************************************
# Load and process income file
muc4a = pd.read_csv(os.path.join(data_path, 'muc4a.csv'))

# Define the required columns (household identifiers + selected variables)
# m4ac6: working months annually (first),  m4ac7: average days work per month (first), m4ac8: average hours work per day (first)
# m4ac16: working months annually (secondary),  m4ac17: average days work per month (secondary), m4ac18: average hours work per day (secondary)
# m4ac21: cash received from secondary job, m4ac22f: other salary 1, m4ac25: other salary 2
# m4ac11: cash received from main job, m4ac12f: other salary
# m4ac21: cash received from secondary job, m4ac22f: other salary 1, m4ac25: other salary 2
columns4a = ['tinh', 'huyen', 'xa', 'diaban', 'hoso', 'matv', 'm4ac6', 'm4ac7', 'm4ac8', 'm4ac11', 
             'm4ac12f', 'm4ac16', 'm4ac17', 'm4ac18', 'm4ac21', 'm4ac22f', 'm4ac25']
muc4a = muc4a[columns4a]

# Create a new column that sums the selected variables for individual's income
muc4a['indi_income'] = muc4a[['m4ac11', 'm4ac12f', 'm4ac21', 'm4ac22f', 'm4ac25']].sum(axis=1)
# Group by household ID to calculate total household income
hh_income = muc4a.groupby(['tinh', 'huyen', 'xa', 'diaban', 'hoso'])['indi_income'].sum().reset_index()
#Rename the column
hh_income.rename(columns={'indi_income': 'HH_Income'}, inplace=True)
# Fill NaN values with zero
muc4a = muc4a.merge(hh_income, on=['tinh', 'huyen', 'xa', 'diaban', 'hoso'], how='left').fillna(0)
# Merge household income back into the individual-level dataset
df = df.merge(muc4a, on=['tinh', 'huyen', 'xa', 'diaban', 'hoso', 'matv'], how='left')

# ******************** CALCULATE TOTAL HOUSEHOLD WEALTH (ASSETS + DURABLE APPLIANCE) *******************************************
# Function to calculate total household wealth
def merge_wealth(df, file_name, columns, value_cols, new_col_name):
    # Load the wealth dataset with selected columns
    data = pd.read_csv(os.path.join(data_path, file_name))[columns]
    # Fill missing values with zero
    data.fillna(0, inplace=True)
    # Calculate household wealth
    if 'm6ac7' in data.columns:  # If percentage ownership column exists (for fixed assets)
        data['m6ac7'] = data['m6ac7'] / 100  # Convert percentage to decimal
        data[new_col_name] = data['m6ac3'] * data['m6ac6'] * data['m6ac7']
    else:  # For durable goods (no ownership percentage)
        data[new_col_name] = data['m6bc3'] * data['m6bc6']

    # Group by household to get total wealth per household
    wealth = data.groupby(['tinh', 'huyen', 'xa', 'diaban', 'hoso'])[new_col_name].sum().reset_index()

    # Merge aggregated wealth data into the main dataframe
    return df.merge(wealth, on=['tinh', 'huyen', 'xa', 'diaban', 'hoso'], how='left')

# Merge wealth data from fixed assets (6A)
# Define the required columns (household identifiers + selected variables)
# m6ac3: quantity of the assets, m6ac6: assets' value at current price, m6ac7: percentage of ownership
df = merge_wealth(df, 'muc6a.csv', ['tinh', 'huyen', 'xa', 'diaban', 'hoso', 'm6ac3', 'm6ac6', 'm6ac7'], ['m6ac3', 'm6ac6', 'm6ac7'], 'HH_wealth1')

# Merge wealth data from durable appliances (6B)
# Define the required columns (household identifiers + selected variables)
# m6bc3: quantity of the durable appliance, m6bc6: durable appliance value at current price
df = merge_wealth(df, 'muc6b.csv', ['tinh', 'huyen', 'xa', 'diaban', 'hoso', 'm6bc3', 'm6bc6'], ['m6bc3', 'm6bc6'], 'HH_wealth2')

# Fill missing values with zero and compute total household wealth
df[['HH_wealth1', 'HH_wealth2']] = df[['HH_wealth1', 'HH_wealth2']].fillna(0)
df['HH_Wealth'] = df['HH_wealth1'] + df['HH_wealth2']

# ******************** CALCULATE TOTAL HOUSEHOLD CONSUMPTION (EXPENDITURES) *******************************************
# Function to process and merge expenditure data
def merge_expenditure(df, file_name, columns, expense_cols, new_col_name):
    data = pd.read_csv(os.path.join(data_path, file_name))[columns] # Load the expenditure datasets with selected columns
    exp = data.groupby(['tinh', 'huyen', 'xa', 'diaban', 'hoso'])[expense_cols].sum().reset_index() # Group the columns by household ID
    exp[new_col_name] = exp[expense_cols].sum(axis=1) # Take the sum of the expenses by households
    exp = exp.drop_duplicates(subset=['tinh', 'huyen', 'xa', 'diaban', 'hoso'], keep='first') # Remove any duplicates, only keep a unique expense values for each household ID
    return df.merge(exp, on=['tinh', 'huyen', 'xa', 'diaban', 'hoso'], how='left') # Merge the datasets to the dataframe df

# Merge all expenditure files

# Load household expenditure of food and drinks during holidays (5A1)
# Define the required columns (household identifiers + selected variables)
# m5a1c2b: expense bought, m5a1c3b: expense self supplied or received
df = merge_expenditure(df, 'muc5a1.csv', ['tinh', 'huyen', 'xa', 'diaban', 'hoso', 'm5a1c2b', 'm5a1c3b'], ['m5a1c2b', 'm5a1c3b'], 'HH_exp1')

# Load household daily expenditure on food and drinks (5A2)
# Define the required columns (household identifiers + selected variables)
# m5a2c6: expense bought, m5a2c10: expense self supplied or received
df = merge_expenditure(df, 'muc5a2.csv', ['tinh', 'huyen', 'xa', 'diaban', 'hoso', 'm5a2c6', 'm5a2c10'], ['m5a2c6', 'm5a2c10'], 'HH_exp2')

# Load household daily expenditure on nonfood and others
# Define the required columns (household identifiers + selected variables)
# m5b1c4: expense recieved, m5b1c5: annual expense
df = merge_expenditure(df, 'muc5b1.csv', ['tinh', 'huyen', 'xa', 'diaban', 'hoso', 'm5b1c4', 'm5b1c5'], ['m5b1c4', 'm5b1c5'], 'HH_exp3')

# Load household annual consumption expenditure
# Define the required columns (household identifiers + selected variables)
# m5b2c2: expense bought, m5b2c3: expense self supplied or received
df = merge_expenditure(df, 'muc5b2.csv', ['tinh', 'huyen', 'xa', 'diaban', 'hoso', 'm5b2c2', 'm5b2c3'], ['m5b2c2', 'm5b2c3'], 'HH_exp4')

# Load other spending that is considered as household expenditure
# Define the required columns (household identifiers + selected variables)
# m5b3c2: annual expense
df = merge_expenditure(df, 'muc5b3.csv', ['tinh', 'huyen', 'xa', 'diaban', 'hoso', 'm5b3c2'], ['m5b3c2'], 'HH_exp5')

# Load household's accomodation expenditure
# Define the required columns (household identifiers + selected variables)
# m7c32: annual water expense, m7c36: annual electricity expense , m7c39: annual garbage collection expense (Expenses)
muc7 = pd.read_csv(os.path.join(data_path, 'muc7.csv'))
columns7 = ['tinh', 'huyen', 'xa', 'diaban', 'hoso', 'm7c32', 'm7c36', 'm7c39']
muc7 = muc7[columns7]
muc7['HH_exp6'] = muc7[['m7c32', 'm7c36', 'm7c39']].sum(axis=1)
muc7.fillna(0, inplace=True)
df = df.merge(muc7, on=['tinh', 'huyen', 'xa', 'diaban', 'hoso'], how='left')

# Aggregate the household income earned from wage  
df['HH_income'] = df['HH_Income'] 

# Aggregate total household consumption expenditure (sum of HH_exp1 to HH_exp6)
df['HH_consumption'] = (
    df['HH_exp1'] + df['HH_exp2'] + df['HH_exp3'] +
    df['HH_exp4'] + df['HH_exp5'] + df['HH_exp6']
)

# Calculate the average household consumption 
df['HH_consumption_avr'] = df['HH_consumption']/ df['hsize']

# Create another new dataframe with only aggregated household income, consumption and wealth
data = df[['tinh', 'huyen', 'xa', 'diaban', 'hoso', 'matv', 'hsize', 'm1ac2',
           'm1ac3', 'm1ac5', 'HH_income', 'HH_consumption', 'HH_consumption_avr', 'HH_Wealth']].copy()

# Rename 'm1ac5' to 'age'
data.rename(columns={'m1ac5': 'age'}, inplace=True)

# Save output to csv file
output_path = r'C:\Users\Do Thu An\OneDrive\Desktop\Dynamic Macroeconomics\Problem sets\Dynamic-Macroeconomics\PS2_Code\Household Model\df_hh.csv'
output_path2 = r'C:\Users\Do Thu An\OneDrive\Desktop\Dynamic Macroeconomics\Problem sets\Dynamic-Macroeconomics\PS2_Code\Household Model\data_processed.csv'
data.to_csv(output_path2, index=False, encoding='utf-8-sig')
df.to_csv(output_path, index=False, encoding='utf-8-sig')