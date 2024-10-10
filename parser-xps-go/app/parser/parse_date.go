package parser

import (
	"fmt"
	"strconv"
	"strings"
	"time"
)

func ParseDate(splittedDate []string) time.Time {

	parsedDay, err := strconv.ParseInt(strings.ReplaceAll(splittedDate[1], ",", ""), 10, 32)
	if err != nil {
		panic(err.Error())
	}

	parsedTime, err := time.Parse("January 2006", fmt.Sprintf("%s %s", splittedDate[0], splittedDate[2]))
	if err != nil {
		panic(err.Error())
	}

	return time.Date(parsedTime.Year(), parsedTime.Month(), int(parsedDay), 12, 0, 0, 0, time.UTC)

}
