package utils

import (
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
