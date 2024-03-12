import pandas as pd
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

# Read the CSV file
data = pd.read_csv('intake-calories.csv')

# Calculate the cumulative total of calories
data['cumulative_calories'] = data['calorie_intake'].cumsum()

# Convert weight to float and handle missing values
data['weight'] = pd.to_numeric(data['weight'], errors='coerce')

# Define the objective function to minimize
def objective_function(params):
    slope, intercept = params
    predicted_weight = (data['cumulative_calories'] - slope * np.arange(len(data))) / 3500 + intercept
    mask = ~data['weight'].isna()
    return np.sum((data.loc[mask, 'weight'] - predicted_weight[mask])**2)

# Initialize the slope (BMR) and intercept
initial_params = [1500, data['weight'].dropna().iloc[0]]

# Minimize the objective function to find the optimal slope (BMR) and intercept
result = minimize(objective_function, initial_params, method='Nelder-Mead')
bmr, initial_weight = result.x

# Calculate the predicted weight using the optimized BMR and intercept
data['predicted_weight'] = (data['cumulative_calories'] - bmr * np.arange(len(data))) / 3500 + initial_weight

# Create a line plot of predicted and actual weight
plt.figure(figsize=(10, 6))
plt.plot(data.index, data['predicted_weight'], label='Predicted Weight')
plt.plot(data.index, data['weight'], 'o', label='Actual Weight')
plt.xlabel('Day')
plt.ylabel('Weight (lbs)')
plt.title('Predicted vs Actual Weight')
plt.legend()
plt.tight_layout()
plt.savefig('weight_prediction.png')
plt.show()

print(f"Estimated BMR: {bmr:.2f} calories per day")
