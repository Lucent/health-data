package models

import "time"

type SingleMeal struct {
	MealName string
	MealType int
	MealDate time.Time
	Foods    SingleFoodProduct
}

type SingleFoodProduct struct {
	FoodName    string
	BrandName   string
	ServingSize string
}
