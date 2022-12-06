package models

import "time"

type SingleMeal struct {
	MealName string
	MealType MealType
	MealDate time.Time
	Food     SingleFoodProduct
}

type SingleFoodProduct struct {
	FoodName    string
	BrandName   string
	ServingSize string
}

// ENUM(Breakfast, Lunch, Dinner, Snacks, Supper)
type MealType int
