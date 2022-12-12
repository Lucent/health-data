package utils

import (
	"bufio"
	"os"
	"strings"

	"golang.org/x/exp/slices"
	"table.reader.lucent/app/constants"
)

func LineCleaner(file *os.File) ([]string, []int, error) {
	tmpStrings := make([]string, 0)
	diaryDateIndex := make([]int, 0)
	counter := 0
	scanner := bufio.NewScanner(file)

	for scanner.Scan() {
		text := scanner.Text() + "\n"
		if strings.Contains(text, "Lucent") || strings.Contains(text, "From") ||
			strings.Contains(text, "Diary") || strings.Contains(text, "diary") || len(scanner.Text()) == 0 {
			continue
		}
		trimmedString := strings.TrimSpace(text)

		splitted := strings.Split(trimmedString, " ")
		if slices.Contains(constants.Months, splitted[0]) && !strings.Contains(splitted[0], "Mayfield") &&
			!strings.Contains(splitted[0], "Octoberfest") {
			diaryDateIndex = append(diaryDateIndex, counter)
		}
		tmpStrings = append(tmpStrings, text)
		counter++
	}
	if err := scanner.Err(); err != nil {
		return []string{}, []int{}, err
	}
	return tmpStrings, diaryDateIndex, nil
}
