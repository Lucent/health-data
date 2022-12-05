package main

import (
	"log"
	"os"

	"table.reader.lucent/app/parser"
	"table.reader.lucent/utils"
)

func main() {
	file, err := os.Open("Food-2012.txt")
	if err != nil {
		log.Fatal(err)
	}
	defer file.Close()
	// client := gosseract.NewClient()
	// defer client.Close()

	// f, err := os.Open("./merged-pics")
	// if err != nil {
	// 	fmt.Println(err)
	// 	return
	// }
	// files, err := f.Readdir(0)
	// if err != nil {
	// 	fmt.Println(err)
	// 	return
	// }

	// var sortedFiles []string

	// for x := range files {
	// 	sortedFiles = append(sortedFiles, files[x].Name())
	// }
	// sort.Slice(sortedFiles, func(i, j int) bool {
	// 	return sortByNumber(sortedFiles[i], sortedFiles[j])
	// })
	// var readFiles []string

	// // Loop through the slice of os.FileInfo objects
	// for y := range sortedFiles {

	// 	// Print the name of each file
	// 	// fmt.Println(file.Name())
	// 	err := client.SetImage("./merged-pics/" + files[y].Name())
	// 	if err != nil {
	// 		panic(err)
	// 		// return
	// 	}
	// 	text, err := client.Text()
	// 	if err != nil {
	// 		panic(err)
	// 	}

	// 	splittedByNewLine := strings.Split(text, "\n")
	// 	for x := range splittedByNewLine {
	// 		readFiles = append(readFiles, splittedByNewLine[x])

	// 	}

	// }

	cleanedFile, dateIndex, err := utils.LineCleaner(file)
	if err != nil {
		log.Fatalln(err.Error())
	}
	for x, v := range dateIndex {
		if x == len(dateIndex)-1 {
			parser.ParseMeal(cleanedFile[v:])
			continue
		}
		currentMeal := cleanedFile[v:dateIndex[x+1]]
		parser.ParseMeal(currentMeal)

	}
}

// // sortByNumber is a function that compares two strings and returns
// // whether the first string should be sorted before or after the second string,
// // based on the numeric value in the strings.
// func sortByNumber(a, b string) bool {
// 	// Get the numeric part of the string by splitting the string on "-" and
// 	// taking the second part of the resulting slice.
// 	numA := strings.ReplaceAll(strings.Split(a, "-")[1], ".png", "")
// 	numB := strings.ReplaceAll(strings.Split(b, "-")[1], ".png", "")

// 	// Convert the numeric string to an int using the `strconv.Atoi` function.
// 	nA, _ := strconv.Atoi(numA)
// 	nB, _ := strconv.Atoi(numB)

// 	// Return whether the first number should be sorted before or after the second number.
// 	return nA < nB
// }
