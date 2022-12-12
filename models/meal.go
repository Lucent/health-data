package models

import "time"

type SingleMeal struct {
	MealName string            `json:"meal_name" bson:"meal_name"`
	MealType MealType          `json:"meal_type" bson:"meal_type"`
	MealDate time.Time         `json:"meal_date" bson:"meal_date"`
	Food     SingleFoodProduct `json:"food" bson:"food"`
}

type SingleFoodProduct struct {
	FoodName    string `json:"food_name" bson:"food_name"`
	BrandName   string `json:"brand_name" bson:"brand_name"`
	ServingSize string `json:"serving_size" bson:"serving_size"`
}

// ENUM(Breakfast, Lunch, Dinner, Snacks, Supper)
type MealType int
