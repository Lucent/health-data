package parser

import (
	"fmt"
	"log"
	"regexp"
	"strings"

	"github.com/samber/lo"
	"table.reader.lucent/app/client"
	"table.reader.lucent/models"
	"table.reader.lucent/utils"
)

func ParseMeal(meal []string) {
	splittedDate := strings.Split(strings.TrimSpace(meal[0]), " ")
	date := ParseDate(splittedDate)
	mealWithoutDate := meal[1:]
	currentMealType := models.MealTypeBreakfast
	counter := 0
	for x := range mealWithoutDate {
		if counter == 5 {
			break
		}
		singleMeal := &models.SingleMeal{}
		var foodItem models.SingleFoodProduct
		singleMeal.MealDate = date
		trimmedString := strings.TrimSpace(mealWithoutDate[x])
		if utils.IsMeal(trimmedString) {
			mealType, err := models.ParseMealType(trimmedString)
			if err != nil {
				log.Fatalln(err)
				return
			}
			currentMealType = mealType
		} else if utils.SkipLine(trimmedString) {
			continue
		}
		lo.TryCatch(func() error {
			foodItem = getFoodDetails(trimmedString)
			return nil
		}, func() {
			fmt.Println("Failed at Line: " + trimmedString)
		})

		if len(foodItem.FoodName) > 0 {
			singleMeal.Food = foodItem
			singleMeal.MealType = currentMealType
			singleMeal.MealName = currentMealType.String()
			err := client.SendHistory(singleMeal)
			if err != nil {
				fmt.Println(err.Error())
			}

		}
		singleMeal = nil
		counter++

	}
}

func getFoodDetails(line string) models.SingleFoodProduct {
	before, after, found := strings.Cut(line, " - ")
	var trimmed string
	var servingSize string
	var singleFoodProduct models.SingleFoodProduct
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
	}

	return singleFoodProduct
}
