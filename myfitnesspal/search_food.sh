#!/bin/bash

urlencode() {
    # jq can be used for urlencode
    jq -rn --arg url "$1" '$url|@uri'
}

# Check if an argument was provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 'food to search for'"
    exit 1
fi

QUERY_STRING=$(urlencode "$1")

curl "https://www.myfitnesspal.com/api/nutrition?query=$QUERY_STRING&page=1" \
  -H 'authority: www.myfitnesspal.com' \
  -H 'accept: application/json, text/plain, */*' \
  -H 'accept-language: en-US,en;q=0.9' \
  -H 'cache-control: no-cache' \
  -H 'pragma: no-cache' \
  -H $'referer: https://www.myfitnesspal.com/food/calorie-chart-nutrition-facts/domino\'s%20handmade' \
  -H 'sec-ch-ua: "Not/A)Brand";v="99", "Microsoft Edge";v="115", "Chromium";v="115"' \
  -H 'sec-ch-ua-mobile: ?0' \
  -H 'sec-ch-ua-platform: "Windows"' \
  -H 'sec-fetch-dest: empty' \
  -H 'sec-fetch-mode: cors' \
  -H 'sec-fetch-site: same-origin' \
  -H 'user-agent: Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/115.0.0.0 Safari/537.36 Edg/115.0.1901.188' \
  --compressed | python -m json.tool
