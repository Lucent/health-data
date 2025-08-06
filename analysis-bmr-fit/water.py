import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize
import sys

NUM_LINES = 1
CALORIES_PER_POUND = 3500
MIN_SEGMENT_LENGTH = 20  # Minimum number of days for each fit line segment

class GlycogenWaterModel:
    def __init__(self, bmr, glycogen_calories_per_kg, max_glycogen_calories, initial_glycogen_calories):
        self.bmr = bmr
        self.glycogen_calories_per_kg = glycogen_calories_per_kg
        self.max_glycogen_calories = max_glycogen_calories
        self.initial_glycogen_calories = initial_glycogen_calories

    def predict(self, data):
        glycogen_calories = np.zeros(len(data))
        glycogen_calories[0] = self.initial_glycogen_calories

        predicted_fat_lean = np.zeros(len(data))
        predicted_fat_lean[0] = data['weight'][0] if not np.isnan(data['weight'][0]) else np.nan

        for i in range(1, len(data)):
            calorie_balance = data['calorie_intake'][i-1] - self.bmr
            if not np.isnan(data['weight'][i]):
                # If weight is available, use it to update fat_lean
                predicted_fat_lean[i] = data['weight'][i] - glycogen_calories[i-1] / self.glycogen_calories_per_kg
            else:
                # If weight is missing, estimate it using calorie intake and BMR
                predicted_fat_lean[i] = predicted_fat_lean[i-1] + calorie_balance / CALORIES_PER_POUND

            glycogen_calories[i] = glycogen_calories[i-1] + calorie_balance
            glycogen_calories[i] = np.clip(glycogen_calories[i], 0, self.max_glycogen_calories)

        glycogen_water_weight = glycogen_calories / self.glycogen_calories_per_kg
        return predicted_fat_lean, glycogen_water_weight

def predict_weight_components(data, bmr, glycogen_calories_per_kg, max_glycogen_calories, initial_glycogen_calories):
    model = GlycogenWaterModel(bmr, glycogen_calories_per_kg, max_glycogen_calories, initial_glycogen_calories)
    predicted_fat_lean, glycogen_water_weight = model.predict(data)

    return pd.DataFrame({
        'date': data['date'],
        'predicted_fat_lean': predicted_fat_lean,
        'predicted_glycogen_water': glycogen_water_weight
    })

def smoothness_score(predicted_fat_lean):
    # Calculate differences between consecutive predicted fat+lean weights
    diffs = np.diff(predicted_fat_lean)
    # Smoothness is measured by the standard deviation of these differences
    return np.std(diffs)

def objective_function(params, data):
    bmr, glycogen_calories_per_kg, max_glycogen_calories, initial_glycogen_calories = params
    predicted = predict_weight_components(data, bmr, glycogen_calories_per_kg, max_glycogen_calories, initial_glycogen_calories)
    return smoothness_score(predicted['predicted_fat_lean'])

# Read CSV from command line argument
input_file = sys.argv[1]
data = pd.read_csv(input_file)

# Convert date column to datetime
data['date'] = pd.to_datetime(data['date'])

# Ensure that no days have missing calorie intake values
if data['calorie_intake'].isna().any():
    raise ValueError("Calorie intake data contains missing values. Please provide complete data.")

# Define initial parameter values and bounds for optimization
initial_params = [2000, 1000, 2000, 1000]  # Initial guesses for [BMR, glycogen_calories_per_kg, max_glycogen_calories, initial_glycogen_calories]
bounds = [(1500, 3000), (500, 1500), (1000, 3000), (500, 2000)]  # Bounds for each parameter

# Optimize parameters to minimize smoothness score
result = minimize(objective_function, initial_params, args=(data,), bounds=bounds, method='L-BFGS-B')
optimized_params = result.x

# Print optimized parameters
print(f"Optimized Parameters:")
print(f"  BMR: {optimized_params[0]:.2f}")
print(f"  Glycogen Calories per kg: {optimized_params[1]:.2f}")
print(f"  Max Glycogen Calories: {optimized_params[2]:.2f}")
print(f"  Initial Glycogen Calories: {optimized_params[3]:.2f}")

# Apply the glycogen and water model with optimized parameters
bmr, glycogen_calories_per_kg, max_glycogen_calories, initial_glycogen_calories = optimized_params
predicted_components = predict_weight_components(data, bmr, glycogen_calories_per_kg, max_glycogen_calories, initial_glycogen_calories)

# Add the predicted components to the original DataFrame
data['predicted_fat_lean'] = predicted_components['predicted_fat_lean']
data['predicted_glycogen_water'] = predicted_components['predicted_glycogen_water']

# Plotting the results
plt.figure(figsize=(10, 6))

# Plot the original weight
plt.plot(data['date'], data['weight'], 'o-', label='Scale Weight', color='blue', alpha=0.8)

# Plot the predicted fat+lean weight
plt.plot(data['date'], data['predicted_fat_lean'], label='Predicted Fat+Lean Weight', color='red', linewidth=2)

# Plot the predicted total weight (fat+lean+glycogen+water)
predicted_total_weight = data['predicted_fat_lean'] + data['predicted_glycogen_water']
plt.plot(data['date'], predicted_total_weight, label='Predicted Total Weight', color='green', linewidth=2)

# Fill the area between predicted fat+lean weight and total weight to show glycogen+water weight
plt.fill_between(data['date'], data['predicted_fat_lean'], predicted_total_weight, color='lightblue', alpha=0.4, label='Predicted Glycogen+Water Weight')

plt.xlabel('Date')
plt.ylabel('Weight (lbs)')
plt.title('Daily Weight with Predicted Fat+Lean and Glycogen+Water Weight')
plt.legend()
plt.grid(True)
plt.show()
