package utils

import (
	"io/fs"
	"log"
	"os"
	"strings"

	"golang.org/x/exp/slices"
	"table.reader.lucent/app/constants"
)

func IsMeal(trimmedString string) bool {
	return slices.Contains(constants.Meals, trimmedString)
}

func SkipLine(trimmedString string) bool {
	return strings.Contains(trimmedString, "TOTAL") ||
		strings.Contains(trimmedString, "Calories Carbs Fat Protein Cholest") ||
		strings.Contains(trimmedString, "EXERCISES") ||
		strings.Contains(trimmedString, "Cardiovascular") ||
		len(trimmedString) < 63
}

func ListFoodHistories() ([]fs.FileInfo, error) {
	f, err := os.Open("./food-history")
	if err != nil {
		log.Println(err)
		return []fs.FileInfo{}, err
	}
	files, err := f.Readdir(0)
	if err != nil {
		log.Println(err)
		return []fs.FileInfo{}, err
	}

	return files, nil
}
