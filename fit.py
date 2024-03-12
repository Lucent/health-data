import pandas as pd
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt

def objective_function(params, data, start_index, end_index):
    slope, intercept = params
    predicted_weight = (data['cumulative_calories'][start_index:end_index+1] - slope * np.arange(end_index-start_index+1)) / 3500 + intercept
    actual_weight = data['weight'][start_index:end_index+1]
    mask = ~actual_weight.isna()
    return np.sum((actual_weight[mask] - predicted_weight[mask])**2)

def find_best_fit_lines(data):
    best_split_day = None
    best_fit_lines = None
    best_total_error = np.inf

    for split_day in range(1, len(data)-1):
        start_data = data[:split_day+1]
        end_data = data[split_day:]

        start_initial_params = [1500, start_data['weight'].dropna().iloc[0] if len(start_data['weight'].dropna()) > 0 else 0]
        end_initial_params = [1500, end_data['weight'].dropna().iloc[0] if len(end_data['weight'].dropna()) > 0 else 0]

        start_result = minimize(objective_function, start_initial_params, args=(data, 0, split_day), method='Nelder-Mead')
        end_result = minimize(objective_function, end_initial_params, args=(data, split_day, len(data)-1), method='Nelder-Mead')

        start_bmr, start_initial_weight = start_result.x
        end_bmr, end_initial_weight = end_result.x

        start_predicted_weight = (start_data['cumulative_calories'] - start_bmr * np.arange(split_day+1)) / 3500 + start_initial_weight
        end_predicted_weight = (end_data['cumulative_calories'] - end_bmr * np.arange(len(data)-split_day)) / 3500 + end_initial_weight

        total_error = start_result.fun + end_result.fun

        if total_error < best_total_error:
            best_split_day = split_day
            best_fit_lines = ((start_bmr, start_initial_weight), (end_bmr, end_initial_weight))
            best_total_error = total_error

    return best_split_day, best_fit_lines

# Read the CSV file
data = pd.read_csv('intake-calories.csv')

# Calculate the cumulative total of calories
data['cumulative_calories'] = data['calorie_intake'].cumsum()

# Convert weight to float and handle missing values
data['weight'] = pd.to_numeric(data['weight'], errors='coerce')

# Find the best fit lines
best_split_day, best_fit_lines = find_best_fit_lines(data)

# Calculate the predicted weights using the best fit lines
(start_bmr, start_initial_weight), (end_bmr, end_initial_weight) = best_fit_lines
start_predicted_weight = (data['cumulative_calories'][:best_split_day+1] - start_bmr * np.arange(best_split_day+1)) / 3500 + start_initial_weight
end_predicted_weight = (data['cumulative_calories'][best_split_day:] - end_bmr * np.arange(len(data)-best_split_day)) / 3500 + end_initial_weight

# Create a line plot of predicted and actual weight
plt.figure(figsize=(10, 6))
plt.plot(data.index[:best_split_day+1], start_predicted_weight, label='Start Fit Line')
plt.plot(data.index[best_split_day:], end_predicted_weight, label='End Fit Line')
plt.plot(data.index, data['weight'], 'o', label='Actual Weight')
plt.xlabel('Day')
plt.ylabel('Weight (lbs)')
plt.title('Best Fit Lines')
plt.legend()
plt.tight_layout()
plt.savefig('best_fit_lines.png')
plt.show()

print(f"Best split day: {best_split_day}")
print(f"Start BMR: {start_bmr:.2f} calories per day")
print(f"End BMR: {end_bmr:.2f} calories per day")
