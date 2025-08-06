# https://www.activtrax.com/workouts.php
# let data = Array.from(document.querySelectorAll("tr[id^='woblock']")).map(row => {
#   let cells = row.querySelectorAll("td");
#   if (cells.length >= 3)
#       return cells[3].innerText.trim();
#   return null;
# }).filter(text => text !== null);
# console.log(JSON.stringify(data));

import json
import os
from datetime import datetime

# Load the data from the JSON file
with open('data.json', 'r') as file:
    data = json.load(file)

missing_pdfs = []

for entry in data:
    try:
        # Extract the relevant part before "at" and split by space to handle date and time separately
        date_str = entry.split('at')[0].strip()
        # Remove the time part to focus on date
        date_str = ' '.join(date_str.split()[:3])
        # Parse the date
        date_obj = datetime.strptime(date_str, "%b %d %Y")
        # Format the date as desired for the PDF filename
        pdf_filename = date_obj.strftime("%Y-%m-%d") + ".pdf"

        if not os.path.isfile(pdf_filename):
            missing_pdfs.append(pdf_filename)
    except Exception as e:
        print(f"Error processing entry '{entry}': {e}")

if missing_pdfs:
    print("Missing PDF files:")
    for pdf in missing_pdfs:
        print(pdf)
else:
    print("All PDF files are present.")

