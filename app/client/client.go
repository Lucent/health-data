package client

import (
	"bytes"
	"encoding/json"
	"errors"
	"fmt"
	"io"
	"log"
	"net/http"

	"table.reader.lucent/models"
)

func doRequest(req *http.Request) (io.ReadCloser, http.Header, error) {
	client := &http.Client{
		// Transport: &http2.Transport{},
	}
	resp, err := client.Do(req)
	if err != nil {
		return nil, http.Header{}, errors.New("request failed")
	}
	fmt.Println(resp.StatusCode)

	return resp.Body, resp.Header, nil
}

func SendHistory(meal *models.SingleMeal) error {
	body, err := json.Marshal(meal)
	if err != nil {
		return err
	}
	req, err := http.NewRequest(http.MethodPost, "http://localhost:4700/api/v1/foods/histories", bytes.NewBuffer(body))
	if err != nil {
		log.Fatal(err)
	}

	_, _, err = doRequest(req)
	if err != nil {
		return err
	}

	return nil
}
