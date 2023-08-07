# Fetching minerals from MyFitnessPal API

https://www.myfitnesspal.com/food/diary/Lucent?date=2023-08-07

This page shows that day's intake. Data sticks around about 4 years. Supposedly it's purged after 2, but we're lucky. Each item is from a database that contains much more, but only macros are displayed. They can be used as a checksum to find the true underlying food data by using the name to query against their database.

`./search_food.sh "Califia Farms - Mocha Iced Coffee with Almondmilk"`

The first JSON result has ID 2321635067 and energy 130 calories for its default serving size of 12 fl oz. I drank 6.5 ounces. 6.5 / 12 * 130 = 70.4, a match! So it's probably safe to use the extra data hidden: iron, potassium, added sugars, saturated fat, sodium, vitamin_d.

It's worth digging deeper and multiplying all macros to see if they match across the board before using an item since this database is littered with near duplicate foods, many with low quality data. I am pretty careful about picking among them.

I'd store the food's full, exact name "Brand - Product" with its JSON entry since it'll likely appear many times.
