import pandas as pd
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import sys

CALORIES_PER_POUND = 3500

def objective_function(params, data, start_index, end_index):
	slope, intercept = params
	predicted_weight = (data['cumulative_calories'][start_index:end_index+1] - slope * np.arange(end_index-start_index+1)) / CALORIES_PER_POUND + intercept
	actual_weight = data['weight'][start_index:end_index+1]
	mask = ~actual_weight.isna()
	return np.sum((actual_weight[mask] - predicted_weight[mask])**2)

def find_best_fit_line(data, start_index, end_index):
	initial_params = [1500, data['weight'][start_index:end_index+1].dropna().iloc[0] if len(data['weight'][start_index:end_index+1].dropna()) > 0 else 0]
	result = minimize(objective_function, initial_params, args=(data, start_index, end_index), method='Nelder-Mead')
	return result.x, result.fun

def find_best_fit_lines(data, num_lines):
	n = len(data)
	dp = np.full((num_lines, n), np.inf)
	path = np.full((num_lines, n), -1, dtype=int)

	for i in range(n):
		line_params, error = find_best_fit_line(data, 0, i)
		dp[0, i] = error

	for j in range(1, num_lines):
		for i in range(j, n):
			for k in range(j-1, i):
				prev_error = dp[j-1, k]
				line_params, current_error = find_best_fit_line(data, k+1, i)
				total_error = prev_error + current_error
				if total_error < dp[j, i]:
					dp[j, i] = total_error
					path[j, i] = k

	split_points = [-1] * num_lines
	split_points[-1] = n-1
	for j in range(num_lines-2, -1, -1):
		split_points[j] = path[j+1, split_points[j+1]]

	return split_points

def plot_best_fit_lines(data, split_points, output_file=None):
	plt.figure(figsize=(10, 6))
	colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
	prev_split = 0

	for i, split_point in enumerate(split_points):
		bmr, initial_weight = find_best_fit_line(data, prev_split, split_point)[0]
		predicted_weight = (data['cumulative_calories'][prev_split:split_point+1] - bmr * np.arange(split_point-prev_split+1)) / CALORIES_PER_POUND + initial_weight
		plt.plot(data.index[prev_split:split_point+1], predicted_weight, color=colors[i%len(colors)], label=f'Fit Line {i+1} (BMR: {bmr:.2f})')
		prev_split = split_point + 1
		print(f"BMR for segment {i+1}: {bmr:.2f}")

	plt.plot(data.index, data['weight'], 'o', label='Actual Weight')

	plt.xlabel('Day')
	plt.ylabel('Weight (lbs)')
	plt.title('Best Fit Lines')
	plt.legend(labelspacing=1.2)
	plt.tight_layout()

	if output_file:
		plt.savefig(output_file)
	else:
		plt.show()

# Read the CSV file
data = pd.read_csv('intake-calories.csv')

# Calculate the cumulative total of calories
data['cumulative_calories'] = data['calorie_intake'].cumsum()

# Convert weight to float and handle missing values
data['weight'] = pd.to_numeric(data['weight'], errors='coerce')

# Specify the number of fit lines
num_lines = 3

# Find the best fit lines
split_points = find_best_fit_lines(data, num_lines)

# Check if an output file is provided as a command line argument
output_file = "best_fit_lines.png"
if len(sys.argv) > 1:
	output_file = sys.argv[1]

# Plot the best fit lines
plot_best_fit_lines(data, split_points, output_file)
