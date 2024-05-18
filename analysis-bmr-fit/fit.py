"""
* Read CSV containing daily calorie intake (NO GAPS) and weight data (intermittent)
* Fill in missing weight values using RMR and linear interpolation
* Split data into segments and finds NUM_LINES-1 breakpoints that gives the lowest total R^2 between new left and right fit lines
"""

import pandas as pd
import numpy as np
from scipy.optimize import minimize
import matplotlib.pyplot as plt
import sys

NUM_LINES = 4
CALORIES_PER_POUND = 3500
MIN_SEGMENT_LENGTH = 20  # Minimum number of days for each fit line segment
OUTPUT_IMAGE = "best_fit_lines.png"

def objective_function(params, data, start_index, end_index):
	slope, intercept = params
	predicted_weight = (data['cumulative_calories'][start_index:end_index+1] - slope * np.arange(end_index-start_index+1)) / CALORIES_PER_POUND + intercept
	actual_weight = data['smoothed_weight'][start_index:end_index+1]
#	actual_weight = data['weight'][start_index:end_index+1]
	return np.sum((actual_weight - predicted_weight)**2)

def fill_missing_weights(data):
	filled_data = data.copy()
	for i in range(len(filled_data)):
		if np.isnan(filled_data['weight'][i]):
			prev_weight_index, next_weight_index = find_nearest_weights(filled_data, i)
			estimated_weight = estimate_weight(filled_data, i, prev_weight_index, next_weight_index)
			filled_data.at[i, 'weight'] = estimated_weight
	return filled_data

def find_nearest_weights(data, index):
	prev_weight_index = next((i for i in range(index - 1, -1, -1) if not np.isnan(data['weight'][i])), None)
	next_weight_index = next((i for i in range(index + 1, len(data)) if not np.isnan(data['weight'][i])), None)
	return prev_weight_index, next_weight_index

def estimate_weight(data, index, prev_weight_index, next_weight_index):
	if prev_weight_index is None or next_weight_index is None:
		return None

	prev_weight = data['weight'][prev_weight_index]
	next_weight = data['weight'][next_weight_index]
	calories_consumed = data['cumulative_calories'][next_weight_index] - data['cumulative_calories'][prev_weight_index]
	time_diff = next_weight_index - prev_weight_index
	weight_diff = (next_weight - prev_weight) * CALORIES_PER_POUND

	BMR = (calories_consumed - weight_diff) / time_diff
	weight_change = (data['cumulative_calories'][index] - data['cumulative_calories'][prev_weight_index]) - BMR * (index - prev_weight_index)
	estimated_weight = prev_weight + weight_change / CALORIES_PER_POUND
	return estimated_weight

def find_best_breakpoint(data, start_index, end_index):
	best_breakpoint = -1
	max_weighted_avg_r2 = -np.inf

	print(f"Look for best breakpoint between {start_index}-{end_index}")
	for i in range(start_index + MIN_SEGMENT_LENGTH - 1, end_index - MIN_SEGMENT_LENGTH + 1):
		left_params, left_r2 = find_best_fit_line(data, start_index, i)
		right_params, right_r2 = find_best_fit_line(data, i, end_index)

		left_segment_length = i - start_index + 1
		right_segment_length = end_index - i
#		weighted_avg_r2 = (left_segment_length * left_r2 + right_segment_length * right_r2) / (left_segment_length + right_segment_length)
		weighted_avg_r2 = left_r2 + right_r2
		print(f"\tcheck {i}: r^2 left {left_r2:.3f}, r^2 right {right_r2:.3f}")

		if weighted_avg_r2 > max_weighted_avg_r2:
			max_weighted_avg_r2 = weighted_avg_r2
			best_breakpoint = i
	
	print(f"\tbest contender: {best_breakpoint}")
	return best_breakpoint, max_weighted_avg_r2

def find_best_fit_line(data, start_index, end_index):
	if end_index - start_index + 1 < MIN_SEGMENT_LENGTH:
		return None, -np.inf

#	initial_weight = data['weight'][start_index]
	initial_weight = data['smoothed_weight'][start_index]
	initial_params = [2500, initial_weight]
	result = minimize(objective_function, initial_params, args=(data, start_index, end_index), method='Nelder-Mead')

	slope, intercept = result.x
	predicted_weight = (data['cumulative_calories'][start_index:end_index+1] - slope * np.arange(end_index-start_index+1)) / CALORIES_PER_POUND + intercept
	#actual_weight = data['weight'][start_index:end_index+1]
	actual_weight = data['smoothed_weight'][start_index:end_index+1]

	sse = np.sum((actual_weight - predicted_weight)**2)
	sst = np.sum((actual_weight - np.mean(actual_weight))**2)
	r_squared = 1 - (sse / sst)

	return (slope, intercept), r_squared

def find_best_fit_lines_recursive(data, num_lines, split_points):
	if num_lines == 1:
		return

	best_breakpoint = -1
	max_weighted_avg_r2 = -np.inf

	for i in range(len(split_points) - 1):
		start_index = split_points[i]
		end_index = split_points[i + 1]

		if end_index - start_index + 1 >= 2 * MIN_SEGMENT_LENGTH:
			breakpoint, weighted_avg_r2 = find_best_breakpoint(data, start_index, end_index)
			if breakpoint != -1 and weighted_avg_r2 > max_weighted_avg_r2:
				max_weighted_avg_r2 = weighted_avg_r2
				best_breakpoint = breakpoint

	insert_position = next(i for i, val in enumerate(split_points) if val > best_breakpoint)
	print(f"Best breakpoint! {best_breakpoint}")
	split_points.insert(insert_position, best_breakpoint)
	find_best_fit_lines_recursive(data, num_lines - 1, split_points)

def find_best_fit_lines(data, num_lines):
	split_points = [0, len(data) - 1]
	find_best_fit_lines_recursive(data, num_lines, split_points)
	return split_points[1:-1]

def smooth_weights_bidirectional(data, REF_RMR, SMOOTHING_FACTOR_RANGE):
	forward_smoothed_data = smooth_weights(data, REF_RMR, SMOOTHING_FACTOR_RANGE, direction='forward')
	backward_smoothed_data = smooth_weights(data, REF_RMR, SMOOTHING_FACTOR_RANGE, direction='backward')

	bidirectional_smoothed_data = data.copy()
	bidirectional_smoothed_data['smoothed_weight'] = (forward_smoothed_data['smoothed_weight'] + backward_smoothed_data['smoothed_weight']) / 2

	return bidirectional_smoothed_data

def smooth_weights(data, REF_RMR, SMOOTHING_FACTOR_RANGE, direction='forward'):
	smoothed_data = data.copy()
	smoothed_data['smoothed_weight'] = smoothed_data['weight']

	last_smoothed_weight = smoothed_data.loc[0, 'smoothed_weight'] if direction == 'forward' else smoothed_data.loc[len(smoothed_data) - 1, 'smoothed_weight']
	cumulative_weight_change = 0

	iter_range = range(1, len(smoothed_data)) if direction == 'forward' else range(len(smoothed_data) - 2, -1, -1)

	for i in iter_range:
		curr_weight = smoothed_data.loc[i, 'weight']
		curr_intake = smoothed_data.loc[i, 'calorie_intake']

		weight_change = (curr_intake - REF_RMR) / CALORIES_PER_POUND
		cumulative_weight_change += weight_change

		if pd.isna(curr_weight):
			smoothed_data.loc[i, 'smoothed_weight'] = last_smoothed_weight + cumulative_weight_change
		else:
			interpolated_weight = last_smoothed_weight + cumulative_weight_change

			# Calculate the adaptive smoothing factor based on the time gap
			time_gap = (smoothed_data.loc[i, 'date'] - smoothed_data.loc[i - 1, 'date']).days if direction == 'forward' else (smoothed_data.loc[i + 1, 'date'] - smoothed_data.loc[i, 'date']).days
			adaptive_smoothing_factor = SMOOTHING_FACTOR_RANGE[0] + (SMOOTHING_FACTOR_RANGE[1] - SMOOTHING_FACTOR_RANGE[0]) * (time_gap - 1) / (MAX_GAP - 1)

			smoothed_data.loc[i, 'smoothed_weight'] = adaptive_smoothing_factor * curr_weight + (1 - adaptive_smoothing_factor) * interpolated_weight
			last_smoothed_weight = smoothed_data.loc[i, 'smoothed_weight']
			cumulative_weight_change = 0

	return smoothed_data

def plot_best_fit_lines(data, split_points, output_file=None):
	plt.figure(figsize=(30, 10))
	colors = ['#FF9999', '#000033', '#0000FF', '#CCCCCC', '#FF00FF', '#FFFF00', '#000000']
	split_points = [0] + split_points + [len(data) - 1]  # Add start and end indices to split_points

	for i in range(len(split_points) - 1):
		start_index = split_points[i]
		end_index = split_points[i + 1]

		(bmr, initial_weight), r_squared = find_best_fit_line(data, start_index, end_index)
		line_length = end_index - start_index + 1
		predicted_weight = (data['cumulative_calories'][start_index:end_index+1] - bmr * np.arange(end_index-start_index+1)) / CALORIES_PER_POUND + initial_weight
		plt.plot(data['date'][start_index:end_index+1], predicted_weight, color=colors[i%len(colors)], label=f'BMR: {bmr:.0f}, R^2: {r_squared:.2f}')

		start_date = data['date'].iloc[start_index].strftime('%Y-%m-%d')
		end_date = data['date'].iloc[end_index].strftime('%Y-%m-%d')
		print(f"BMR for segment {i+1} ({start_date} to {end_date}): {bmr:.0f}, R^2: {r_squared:.3f}")

	plt.plot(data['date'], data['weight'], '-o', label='Actual Weight')
	plt.plot(data['date'], data['smoothed_weight'], '-', label='Smoothed Weight')

	plt.xlabel('Date')
	plt.ylabel('Weight (lbs)')
	plt.title('Best Fit Lines')
	plt.legend(labelspacing=1.2)
	plt.xticks(rotation=45)
	plt.tight_layout()

	plt.savefig(output_file)


data = pd.read_csv('intake-weight.csv')

data['date'] = pd.to_datetime(data['date'])

data['cumulative_calories'] = data['calorie_intake'].cumsum()

data['weight'] = pd.to_numeric(data['weight'], errors='coerce')

REF_RMR = 2400
SMOOTHING_FACTOR_RANGE = (0.2, 0.8)  # Min and max smoothing factors based on time gap
MAX_GAP = 10  # Maximum expected time gap between weight measurements

smoothed_data = smooth_weights_bidirectional(data, REF_RMR, SMOOTHING_FACTOR_RANGE)
data['smoothed_weight'] = smoothed_data['smoothed_weight']

split_points = find_best_fit_lines(smoothed_data, NUM_LINES)

if len(sys.argv) > 1:
	output_file = sys.argv[1]

plot_best_fit_lines(data, split_points, OUTPUT_IMAGE)
