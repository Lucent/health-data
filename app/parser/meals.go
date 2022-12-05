package parser

import (
	"fmt"
	"regexp"
	"strings"

	"github.com/samber/lo"
	"table.reader.lucent/models"
)

func ParseMeal(meal []string) {
	splittedDate := strings.Split(strings.TrimSpace(meal[0]), " ")

	date := ParseDate(splittedDate)
	fmt.Println("CURRENTLY AT DATE: " + date.String())
	mealWithoutDate := meal[1:]

	for x := range mealWithoutDate {
		trimmedString := strings.TrimSpace(mealWithoutDate[x])
		fmt.Println("orig: " + trimmedString)

		if strings.Contains(trimmedString, "Breakfast") || trimmedString == "Lunch" ||
			trimmedString == "Dinner" ||
			trimmedString == "Snacks" || strings.Contains(trimmedString, "Supper") ||
			strings.Contains(trimmedString, "TOTAL") ||
			trimmedString == "FOODS                                                                          Calories Carbs Fat Protein Cholest      Sodium Sugars Fiber" ||
			strings.Contains(trimmedString, "EXERCISES") ||
			strings.Contains(trimmedString, "Cardiovascular") ||

			len(trimmedString) < 63 {
			fmt.Println("LOST :" + trimmedString)
			continue
		}

		foodItem := getFoodDetails(trimmedString)
		fmt.Println("brand: "+foodItem.BrandName, "food: "+foodItem.FoodName)
	}
}

func getFoodDetails(line string) models.SingleFoodProduct {
	before, after, found := strings.Cut(line, " - ")
	var trimmed string
	var servingSize string
	var singleFoodProduct models.SingleFoodProduct
	lo.TryCatch(func() error {
		if found {
			lastIndex := strings.Index(after, ",")
			space := regexp.MustCompile(`\s+`)
			digits := regexp.MustCompile(`[^(][0-9]{3}[^cup][^oz]`)

			if lastIndex != -1 {
				trimmed = strings.TrimSpace(after[:lastIndex])
				s := space.ReplaceAllString(after, " ")
				result := digits.Split(s, -1)
				splitted := strings.Split(result[0], ",")
				digits = regexp.MustCompile(`[^(][0-9]{2}[^cup][^oz][^Pizza][^Tbsp]`)
				servingSize = strings.TrimSpace(digits.Split(splitted[1], 2)[0])
			}
			singleFoodProduct = models.SingleFoodProduct{
				BrandName:   strings.TrimSpace(before),
				FoodName:    trimmed,
				ServingSize: servingSize,
			}
			// return strings.TrimSpace(splitted[0])
		}
		return nil
	}, func() {
		// caught = true

	})

	// if len(parts) == 2 { // Get the string after the dash
	// 	afterDash := strings.TrimSpace(parts[1])

	// 	// Find the index of the first comma in this string
	// 	commaIndex := strings.Index(afterDash, ",")

	// 	// The foodname is the substring between the dash and the comma
	// 	return strings.TrimSpace(afterDash[:commaIndex])
	// } else if len(parts) == 3 {
	// 	afterDash := strings.TrimSpace(parts[1] + " " + parts[2])

	// 	// Find the index of the first comma in this string
	// 	commaIndex := strings.Index(afterDash, ",")

	// 	return strings.TrimSpace(afterDash[:commaIndex])

	// }
	// return ""
	return singleFoodProduct
}
