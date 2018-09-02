import axios from 'axios'
import * as FilePond from 'filepond'
import { saveAs } from 'file-saver/FileSaver'
import { countBy, uniq, zip } from 'lodash'
import csv from 'papaparse'
import pako from 'pako'

$('#mainTab a').on('click', function(e) {
  e.preventDefault()
  $(this).tab('show')
})

$(function() {
  $('[data-toggle="tooltip"]').tooltip()
})

const resultLink = document.getElementById('link-results')

const submitButton = document.getElementById('btn-submit')
submitButton.addEventListener('click', function() {
  resultLink.click()
  run()
})

let data, exampleData, readGroups, summary
const exampleButton = document.getElementById('btn-example')
exampleButton.addEventListener('click', showExample)
const urlExample = document.getElementById('link-example')

const inputFile = document.getElementById('inputFile')
const chartsContainer = document.getElementById('charts-container')
const resultContainer = document.getElementById('result-container')
const resultInfo = document.getElementById('result-info')
const resultError = document.getElementById('result-error')
const selectSample = document.getElementById('select-sample')
const selectReadGroup = document.getElementById('select-rg')
const selectToc = document.getElementById('select-toc')
const summaryTab = document.getElementById('summary-tab')

const fileUpload = FilePond.create(inputFile)

function run() {
  const fileObjects = fileUpload.getFiles()

  if (fileObjects.length === 0) {
    showError('Error: no files specified')
    return
  }

  const fileCounts = countBy(fileObjects.map(f => f.filename))
  for (const [fileName, count] of Object.entries(fileCounts)) {
    if (count > 1) {
      showError(`Error: file <code>${fileName}</code> specified multiple times`)
      return
    }
  }

  hideElement(resultContainer)
  hideElement(resultError)
  showElement(resultInfo)
  summaryTab.innerHTML = ''

  mergeInputs(fileObjects).then(() => {
    const sampleCounts = countBy(data.samples.map(s => s.id))
    for (const [sampleId, count] of Object.entries(sampleCounts)) {
      if (count > 1) {
        showError(`Error: sample <code>${sampleId}</code> is not unique`)
        return
      }
    }
    handleSuccess(data)
  })
}

function readFile(fileObject) {
  const fileReader = new FileReader()
  const isGzip = fileObject.fileExtension === 'gz'

  if (isGzip) {
    fileReader.readAsArrayBuffer(fileObject.file)
  } else {
    fileReader.readAsText(fileObject.file)
  }

  return new Promise(resolve => {
    fileReader.onload = event => {
      let content = event.target.result
      if (isGzip) {
        content = pako.ungzip(content, { to: 'string' })
      }
      data = JSON.parse(content)
      resolve(data)
    }
  })
}

function mergeInputs(fileObjects) {
  const fileReads = []
  for (const fileObject of fileObjects) {
    fileReads.push(readFile(fileObject))
  }
  return Promise.all(fileReads).then(fileData => {
    data = {
      samples: fileData
        .map(d => d.samples)
        .reduce((acc, cur) => acc.concat(cur))
        .sort((s1, s2) => s1.id.localeCompare(s2.id))
    }
    consolidateSummaries(data)
  })
}

function consolidateSummaries(data) {
  const allColumns = uniq(
    data.samples
      .map(s => s.summary.data.columns)
      .reduce((acc, cur) => acc.concat(cur))
  )
  for (const sample of data.samples) {
    const oldColumns = sample.summary.data.columns
    const oldRows = sample.summary.data.rows.map(
      row => new Map(zip(oldColumns, row))
    )
    sample.summary.data.columns = allColumns
    sample.summary.data.rows = []
    for (const row of oldRows) {
      sample.summary.data.rows.push(
        allColumns.map(col => (row.has(col) ? row.get(col) : null))
      )
    }
  }
}

function handleSuccess(data) {
  hideElement(resultInfo)
  hideElement(resultError)
  showElement(resultContainer)

  chartsContainer.innerHTML = ''

  const samples = data.samples.map(sample => sample.id)
  readGroups = {}
  data.samples.forEach(sample => {
    readGroups[sample.id] = sample.readGroups.map(rg => rg.id)
  })

  selectSample.innerHTML = samples.map(s => `<option>${s}</option>`).join('')

  selectReadGroup.innerHTML = readGroups[samples[0]]
    .map(rg => `<option>${rg}</option>`)
    .join('')

  summary = {
    title: data.samples[0].summary.title,
    data: {
      columns: data.samples[0].summary.data.columns,
      rows: data.samples
        .map(s => s.summary.data.rows)
        .reduce((acc, cur) => acc.concat(cur))
    }
  }

  summaryTable(summary, true)
  populateToc(samples[0], readGroups[samples[0]][0])
  vis(data, samples[0], readGroups[samples[0]][0])
}

window.handleReadGroupSelectChange = handleReadGroupSelectChange
function handleReadGroupSelectChange() {
  const sample = selectSample.value
  const readGroup = selectReadGroup.value
  chartsContainer.innerHTML = ''
  populateToc(sample, readGroup)
  vis(data, sample, readGroup)
}

window.handleSampleSelectChange = handleSampleSelectChange
function handleSampleSelectChange() {
  const sample = selectSample.value
  selectReadGroup.innerHTML = readGroups[sample]
    .map(rg => `<option>${rg}</option>`)
    .join('')
  const readGroup = readGroups[sample][0]
  chartsContainer.innerHTML = ''
  populateToc(sample, readGroup)
  vis(data, sample, readGroup)
}

function populateToc(sample, readGroup) {
  const dataRg = data.samples
    .find(s => s.id === sample)
    .readGroups.find(rg => rg.id === readGroup)

  selectToc.innerHTML = `${dataRg.metrics.map(
    metric => `<option value="${metric.id}">${metric.title}</option>`
  )}`
}

window.handleTocChange = handleTocChange
function handleTocChange() {
  const id = selectToc.value
  document.getElementById(id).scrollIntoView()
}

const chartDispatch = {
  bar: chart,
  line: chart,
  table: table
}

function vis(data, sample, readGroup) {
  const dataRg = data.samples
    .find(s => s.id === sample)
    .readGroups.find(rg => rg.id === readGroup)

  for (const metric of dataRg.metrics) {
    chartDispatch[metric.type](metric, chartsContainer)
  }
}

function chart(metricData, parent) {
  const container = document.createElement('div')
  container.id = metricData.id
  parent.appendChild(container)

  const xData = metricData.x.data[0].values
  const chartData = []
  for (const y of metricData.y.data) {
    const trace = {
      x: xData,
      y: y.values,
      name: y.title || ''
    }
    if (metricData.type === 'bar') {
      trace.type = 'bar'
    }
    chartData.push(trace)
  }

  const layout = {
    title: metricData.title,
    xaxis: {
      title: metricData.x.axis.title,
      zeroline: false
    },
    yaxis: {
      title: metricData.y.axis.title,
      zeroline: false
    }
  }

  if (metricData.type === 'bar') {
    layout.barmode = metricData.options.layout
  }

  if (metricData.x.axis.range) {
    layout.xaxis.range = metricData.x.axis.range
  }

  if (metricData.y.axis.range) {
    layout.yaxis.range = metricData.y.axis.range
  }

  Plotly.newPlot(container, chartData, layout)
}

// TODO consolidate / refactor table functions

function table(tableData, parent) {
  const html = `
    <h4>${tableData.title}</h4>
    <div style="overflow-x: auto;">
      <table class="table table-sm table-striped table-hover">
        <thead>
          <tr>
            ${tableData.data.columns
              .map(title => `<th scope="col">${title}</th>`)
              .join('')}
          </tr>
        </thead>
        <tbody>
          ${tableData.data.rows
            .map(
              row => `<tr>
              ${row
                .map(
                  (value, i) =>
                    `<td title="${tableData.data.columns[i]}">${value}</td>`
                )
                .join('')}
            </tr>`
            )
            .join('')}
        </tbody>
      </table>
    </div>
  `
  const element = document.createElement('div')
  element.id = tableData.id
  element.innerHTML = html
  parent.appendChild(element)
}

function summaryTable(tableData, transpose = false) {
  let html = `
    <h4>${tableData.title}</h4>
    <div class="mb-2">
      <button type="button" class="btn btn-outline-primary" onclick="transpose()">
        <i class="fas fa-redo-alt" style="margin-right: 5px;"></i>
        Transpose table
      </button>
      <button type="button" class="btn btn-outline-primary" onclick="summaryDownload()">
        <i class="fas fa-file-download" style="margin-right: 5px;"></i>
        Download .csv
      </button>
    </div>
    <div style="overflow-x: auto;">
    `
  if (transpose) {
    const rows = zip(tableData.data.columns, ...tableData.data.rows)
    html += `
      <table id="summary-table" class="table table-sm table-striped table-hover" data-transposed>
        <tbody>
        ${rows
          .map(
            row => `<tr>
            ${row
              .map((value, i) => {
                if (i === 0) {
                  return `<th scope="row">${value}</th>`
                }
                return `<td title="${row[0]}">${value === null ? '—' : value}</td>`
              })
              .join('')}
          </tr>`
          )
          .join('')}
        </tbody>
      </table>
    `
  } else {
    html += `
      <table id="summary-table" class="table table-sm table-striped table-hover">
        <thead>
          <tr>
            ${tableData.data.columns
              .map(title => `<th scope="col">${title}</th>`)
              .join('')}
          </tr>
        </thead>
        <tbody>
          ${tableData.data.rows
            .map(
              row => `<tr>
              ${row
                .map(
                  (value, i) =>
                    `<td title="${tableData.data.columns[i]}">${value === null ? '—' : value}</td>`
                )
                .join('')}
            </tr>`
            )
            .join('')}
        </tbody>
      </table>
    `
  }
  html += '</div>'
  summaryTab.innerHTML = html
}

window.transpose = transpose
function transpose() {
  const tableElement = document.querySelector('#summary-table')
  const isTransposed = 'transposed' in tableElement.dataset
  summaryTable(summary, !isTransposed)
}

window.summaryDownload = summaryDownload
function summaryDownload() {
  const tableElement = document.querySelector('#summary-table')
  const isTransposed = 'transposed' in tableElement.dataset
  let data
  if (isTransposed) {
    data = zip(summary.data.columns, ...summary.data.rows)
  } else {
    data = [summary.data.columns, ...summary.data.rows]
  }
  const csvData = csv.unparse(data)
  const blob = new Blob([csvData], { type: 'text/plain;charset=utf-8' })
  // TODO generate file name from input
  saveAs(blob, 'alfred-summary-stats.csv')
}

function showExample() {
  hideElement(resultContainer)
  hideElement(resultError)
  chartsContainer.innerHTML = ''
  summaryTab.innerHTML = ''
  showElement(resultInfo)
  resultLink.click()

  if (exampleData) {
    data = exampleData
    handleSuccess(data)
  } else {
    axios
      .get(urlExample, {
        responseType: 'arraybuffer'
      })
      .then(response => {
        const content = pako.ungzip(response.data, { to: 'string' })
        exampleData = JSON.parse(content)
        data = exampleData
        handleSuccess(data)
      })
      .catch(error => {
        // FIXME proper error handling
        console.error(error)
      })
  }
}

function showError(message) {
  hideElement(resultContainer)
  hideElement(resultInfo)
  summaryTab.innerHTML = ''

  showElement(resultError)
  document.querySelector('#error-message').innerHTML = message
}

function showElement(element) {
  element.classList.remove('d-none')
}

function hideElement(element) {
  element.classList.add('d-none')
}
