import datetime
import myfitnesspal
import csv

client = myfitnesspal.Client()

start_date = datetime.date(2018, 2, 19)
end_date = datetime.date(2018, 2, 20)
current_date = start_date

csv_file = open('calories.csv', 'a')
csv_writer = csv.writer(csv_file)

csv_writer.writerow(['Date', 'Calories'])

while current_date <= end_date:
	day = client.get_date(current_date.year, current_date.month, current_date.day)
	print(day)
	calories = day.totals['calories']
	csv_writer.writerow([current_date, calories])
	current_date += datetime.timedelta(days=1)

csv_file.close()
