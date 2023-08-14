const fs = require('fs');
const jsdom = require('jsdom');
const mhtml2html = require('mhtml2html');
const { createObjectCsvWriter } = require('csv-writer');
const { JSDOM } = jsdom;

const directoryPath = '../intake';

function convertDateFormat(dateStr) {
	let date = new Date(Date.parse(dateStr));

	let year = date.getFullYear();
	let month = date.getMonth() + 1;  // getMonth() returns month index starting from 0
	let day = date.getDate();

	return `${year}-${month.toString().padStart(2, '0')}-${day.toString().padStart(2, '0')}`;
}

function extractHTML(filePath) {
	const htmlContent = fs.readFileSync(filePath, 'utf-8');
	const dom = new JSDOM(htmlContent);
	const document = dom.window.document;
	return document;
}

function extractMHTML(filePath) {
	const htmlContent = fs.readFileSync(filePath, 'utf-8');
	const dom = mhtml2html.convert(htmlContent, { parseDOM: (html) => new JSDOM(html) });
	const document = dom.window.document;
	return document;
}

function writeDataToCsv(data, outputFilename) {
	const csvWriter = createObjectCsvWriter({
		path: outputFilename,
		header: [
			{ id: 'date', title: 'Date' },
			{ id: 'calories', title: 'Calories' }
		],
		alwaysQuote: false
	});

	csvWriter.writeRecords(data);
}

fs.readdir(directoryPath, { recursive: true }, (err, files) => {
	if (err) {
		console.error('Error reading directory:', err);
		return;
	}

	let pairs = [];
	files.forEach(file => {
		const filePath = `${directoryPath}/${file}`;
		let content;
		if (file.endsWith('.mhtml'))
			content = extractMHTML(filePath);
		else if (file.endsWith('.html'))
			content = extractHTML(filePath);
		else return;
		console.log(filePath);

		if (content.querySelector("main")) { // new
			content.querySelectorAll('main > div:not(:first-of-type)').forEach(day => {
				pairs.push({
					date: convertDateFormat(day.querySelector(":scope > p").textContent),
					calories: day.querySelector(":scope table tbody:last-of-type th:nth-of-type(2)").textContent
				});
			});
		} else { // old
			content.querySelectorAll('h2').forEach(day => {
				pairs.push({
					date: convertDateFormat(day.textContent),
					calories: day.nextElementSibling.querySelector("tfoot td:nth-of-type(2)").textContent.replace(",","")
				});
			});
		}
	});

	writeDataToCsv(pairs, 'output.csv');
	console.log('Data written to output.csv');
});
