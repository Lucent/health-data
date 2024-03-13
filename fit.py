import pandas as pd
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import sys

CALORIES_PER_POUND = 3500
MIN_SEGMENT_LENGTH = 12  # Minimum number of days for each fit line segment

def objective_function(params, data, start_index, end_index):
	slope, intercept = params
	predicted_weight = (data['cumulative_calories'][start_index:end_index+1] - slope * np.arange(end_index-start_index+1)) / CALORIES_PER_POUND + intercept
	actual_weight = data['weight'][start_index:end_index+1]
	mask = ~actual_weight.isna()
	return np.sum((actual_weight[mask] - predicted_weight[mask])**2)

def find_best_fit_line(data, start_index, end_index):
	if end_index - start_index + 1 < MIN_SEGMENT_LENGTH:
		return None, np.inf
	initial_params = [2500, data['weight'][start_index:end_index+1].dropna().iloc[0] if len(data['weight'][start_index:end_index+1].dropna()) > 0 else 0]
	result = minimize(objective_function, initial_params, args=(data, start_index, end_index), method='Nelder-Mead')
	return result.x, result.fun

def find_best_breakpoint(data, start_index, end_index):
	best_breakpoint = -1
	min_total_error = np.inf

	print(f"Look for best breakpoint between {start_index}-{end_index}")
	for i in range(start_index + MIN_SEGMENT_LENGTH - 1, end_index - MIN_SEGMENT_LENGTH + 1):
		left_params, left_error = find_best_fit_line(data, start_index, i)
		right_params, right_error = find_best_fit_line(data, i + 1, end_index)
#		print(f"\tbreakpoint {i} error: {left_error+right_error}")

		if left_params is not None and right_params is not None:
			total_error = left_error + right_error
			if total_error < min_total_error:
				min_total_error = total_error
				best_breakpoint = i

	print(f"\tfound at {best_breakpoint} with error {min_total_error}")
	return best_breakpoint, min_total_error

def find_best_fit_lines_recursive(data, num_lines, split_points):
	if num_lines == 1:
		return

	best_breakpoint = -1
	min_avg_error = np.inf
	best_segment_index = -1

	for i in range(len(split_points) - 1):
		start_index = split_points[i] + 1
		end_index = split_points[i + 1]

		if end_index - start_index + 1 >= 2 * MIN_SEGMENT_LENGTH:
			breakpoint, avg_error = find_best_breakpoint(data, start_index, end_index)
			if breakpoint != -1 and avg_error < min_avg_error:
				min_avg_error = avg_error
				best_breakpoint = breakpoint
				best_segment_index = i

	if best_breakpoint != -1:
		split_points.insert(best_segment_index + 1, best_breakpoint)
		print(f"\twinner {best_breakpoint} had lower error")
		find_best_fit_lines_recursive(data, num_lines - 1, split_points)
	else:
		print("No suitable breakpoint found in any segment")

def find_best_fit_lines(data, num_lines):
	split_points = [0, len(data) - 1]  # Initialize split_points with start and end indices
	find_best_fit_lines_recursive(data, num_lines, split_points)
	return split_points[1:-1]  # Exclude the start and end indices

def plot_best_fit_lines(data, split_points, output_file=None):
	plt.figure(figsize=(10, 6))
	colors = ['r', 'g', 'b', 'c', 'm', 'y', 'k']
	split_points = [0] + split_points + [len(data) - 1]  # Add start and end indices to split_points

	for i in range(len(split_points) - 1):
		start_index = split_points[i]
		end_index = split_points[i + 1]

		(bmr, initial_weight), total_error = find_best_fit_line(data, start_index, end_index)
		line_length = end_index - start_index + 1
		avg_error = total_error / line_length
		predicted_weight = (data['cumulative_calories'][start_index:end_index+1] - bmr * np.arange(end_index-start_index+1)) / CALORIES_PER_POUND + initial_weight
		plt.plot(data['date'][start_index:end_index+1], predicted_weight, color=colors[i%len(colors)], label=f'BMR: {bmr:.0f}, Avg err: {avg_error:.2f}')

		start_date = data['date'].iloc[start_index].strftime('%Y-%m-%d')
		end_date = data['date'].iloc[end_index].strftime('%Y-%m-%d')
		print(f"BMR for segment {i+1} ({start_date} to {end_date}): {bmr:.0f}, Avg Error: {avg_error:.3f}")

	plt.plot(data['date'], data['weight'], 'o', label='Actual Weight')

	plt.xlabel('Date')
	plt.ylabel('Weight (lbs)')
	plt.title('Best Fit Lines')
	plt.legend(labelspacing=1.2)
	plt.xticks(rotation=45)
	plt.tight_layout()

	if output_file:
		plt.savefig(output_file)
	else:
		plt.show()

# Read the CSV file
data = pd.read_csv('intake-weight.csv')

# Convert 'date' column to datetime
data['date'] = pd.to_datetime(data['date'])

# Calculate the cumulative total of calories
data['cumulative_calories'] = data['calorie_intake'].cumsum()

# Convert weight to float and handle missing values
data['weight'] = pd.to_numeric(data['weight'], errors='coerce')

# Specify the number of fit lines
NUM_LINES = 4

# Find the best fit lines
split_points = find_best_fit_lines(data, NUM_LINES)

# Check if an output file is provided as a command line argument
output_file = "best_fit_lines.png"
if len(sys.argv) > 1:
	output_file = sys.argv[1]

# Plot the best fit lines
plot_best_fit_lines(data, split_points, output_file)
