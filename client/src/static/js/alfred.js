import axios from 'axios'
import Choices from 'choices.js'
import * as FilePond from 'filepond'
import { saveAs } from 'file-saver/FileSaver'
import { countBy, uniq, zip } from 'lodash'
import csv from 'papaparse'
import pako from 'pako'

import examples from '../examples/examples.json'

document.addEventListener('DOMContentLoaded', () => {
  // DOM elements
  const summaryTab = document.getElementById('summary-tab')
  const resultLink = document.getElementById('link-results')
  const submitButton = document.getElementById('btn-submit')
  const inputFile = document.getElementById('inputFile')
  const chartsContainer = document.getElementById('charts-container')
  const resultContainer = document.getElementById('result-container')
  const resultInfo = document.getElementById('result-info')
  const resultError = document.getElementById('result-error')
  const selectSample = Select(document.getElementById('select-sample'))
  const selectReadGroup = Select(document.getElementById('select-rg'))
  const selectExamples = Select(document.getElementById('select-examples'))
  const selectToc = Select(document.getElementById('select-toc'))

  // Defensive: log missing elements
  [
    ['resultLink', resultLink],
    ['submitButton', submitButton],
    ['inputFile', inputFile],
    ['chartsContainer', chartsContainer],
    ['resultContainer', resultContainer],
    ['resultInfo', resultInfo],
    ['resultError', resultError],
    ['summaryTab', summaryTab]
  ].forEach(([name, el]) => {
    if (!el) {
      console.warn(`Element '${name}' is missing from DOM!`);
    }
  });

  // Only continue if all main elements exist
  if (!submitButton || !resultLink || !inputFile || !chartsContainer || !resultContainer || !resultInfo || !resultError || !summaryTab) {
    console.error('Missing one or more required elements. Aborting initialization.');
    return;
  }

  $('#mainTab a').on('click', function(e) {
    e.preventDefault()
    $(this).tab('show')
  })

  $('[data-toggle="tooltip"]').tooltip()

  submitButton.addEventListener('click', function() {
    resultLink.click()
    run()
  })

  let data, exampleData, readGroups, summary

  const fileUpload = FilePond.create(inputFile)

  selectExamples.setChoices(
    examples.map(ex => ({
      value: ex.filename,
      label: ex.title
    })),
    'value',
    'label',
    true
  )

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
      if (data) {
        handleSuccess(data)
      }
    })
  }

  function Select(element, options = {}) {
    if (!element) return { setChoices: () => {}, getValue: () => {} }
    return new Choices(
      element,
      Object.assign(
        {
          shouldSort: false,
          searchResultLimit: 50,
          searchFields: ['label']
        },
        options
      )
    )
  }

  function readFile(file) {
    const fileReader = new FileReader()
    const isGzip = file.name.endsWith('.gz')

    if (isGzip) {
      fileReader.readAsArrayBuffer(file)
    } else {
      fileReader.readAsText(file)
    }

    return new Promise((resolve, reject) => {
      fileReader.onload = event => {
        let content = event.target.result
        if (isGzip) {
          content = pako.ungzip(content, { to: 'string' })
        }
        try {
          data = JSON.parse(content)
          if (!data.samples) {
            reject(`Error(${file.name}): wrong format, missing 'samples' key.`)
          }
          resolve(data)
        } catch (error) {
          reject(`Error(${file.name}): not a JSON file.`)
        }
      }
    })
  }

  function mergeInputs(fileObjects) {
    const fileReads = []
    for (const fileObject of fileObjects) {
      fileReads.push(readFile(fileObject.file))
    }
    return Promise.all(fileReads)
      .then(fileData => {
        data = {
          samples: fileData
            .map(d => d.samples)
            .reduce((acc, cur) => acc.concat(cur))
        }
        consolidateSummaries(data)
      })
      .catch(error => {
        data = undefined
        showError(error)
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
    chartsContainer.innerHTML = ''

    const samples = uniq(data.samples.map(sample => sample.id))

    readGroups = {}
    data.samples.forEach(sample => {
      if (!(sample.id in readGroups)) {
        readGroups[sample.id] = new Set()
      }
      for (const rg of sample.readGroups) {
        if (readGroups[sample.id].has(rg.id)) {
          showError(`Error: read groups of sample '${sample.id}' are not unique.`)
        }
        readGroups[sample.id].add(rg.id)
      }
    })

    selectSample.setChoices(
      samples.map((s, i) => ({
        value: s,
        label: s,
        selected: i === 0
      })),
      'value',
      'label',
      true
    )

    selectReadGroup.setChoices(
      [...readGroups[samples[0]].values()].map((rg, i) => ({
        value: rg,
        label: rg,
        selected: i === 0
      })),
      'value',
      'label',
      true
    )

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
    populateToc(samples[0], [...readGroups[samples[0]].values()][0])
    vis(data, samples[0], [...readGroups[samples[0]].values()][0])
  }

  window.handleExampleSelectChange = handleExampleSelectChange
  function handleExampleSelectChange() {
    const example = selectExamples.getValue(true)
    showExample(example)
  }

  window.handleReadGroupSelectChange = handleReadGroupSelectChange
  function handleReadGroupSelectChange() {
    const sample = selectSample.getValue(true)
    const readGroup = selectReadGroup.getValue(true)
    chartsContainer.innerHTML = ''
    populateToc(sample, readGroup)
    showElement(resultInfo)
    setTimeout(() => {
      vis(data, sample, readGroup)
    })
  }

  window.handleSampleSelectChange = handleSampleSelectChange
  function handleSampleSelectChange() {
    const sample = selectSample.getValue(true)
    const rgs = [...readGroups[sample].values()]
    selectReadGroup.setChoices(
      rgs.map((rg, i) => ({
        value: rg,
        label: rg,
        selected: i === 0
      })),
      'value',
      'label',
      true
    )
    const readGroup = rgs[0]
    chartsContainer.innerHTML = ''
    populateToc(sample, readGroup)
    showElement(resultInfo)
    setTimeout(() => {
      vis(data, sample, readGroup)
    })
  }

  function populateToc(sample, readGroup) {
    const dataRg = data.samples
      .filter(s => s.id === sample)
      .find(s => s.readGroups.find(rg => rg.id === readGroup))
      .readGroups.find(rg => rg.id === readGroup)

    selectToc.setChoices(
      dataRg.metrics.map((metric, i) => ({
        value: metric.id,
        label: metric.title,
        selected: i === 0
      })),
      'value',
      'label',
      true
    )
  }

  window.handleTocChange = handleTocChange
  function handleTocChange() {
    const targetId = selectToc.getValue(true)
    setTimeout(() => {
      const el = document.getElementById(targetId)
      if (el) el.scrollIntoView()
    }, 0)
  }

  const chartDispatch = {
    bar: chart,
    line: chart,
    table: table
  }

  function vis(data, sample, readGroup) {
    hideElement(resultInfo)
    hideElement(resultError)
    showElement(resultContainer)

    const dataRg = data.samples
      .filter(s => s.id === sample)
      .find(s => s.readGroups.find(rg => rg.id === readGroup))
      .readGroups.find(rg => rg.id === readGroup)

    for (const metric of dataRg.metrics) {
      chartDispatch[metric.type](metric, chartsContainer)
    }
  }

  function chart(metricData, parent) {
    const container = document.createElement('div')
    container.id = metricData.id
    parent.appendChild(container)

    const chartData = []
    for (let i = 0; i < metricData.y.data.length; i += 1) {
      const xData =
        metricData.x.data.length === metricData.y.data.length
          ? metricData.x.data[i].values
          : metricData.x.data[0].values
      const trace = {
        x: xData,
        y: metricData.y.data[i].values,
        name: metricData.y.data[i].title || ''
      }
      if (metricData.type === 'bar') {
        trace.type = 'bar'
      }
      chartData.push(trace)
    }

    let title = metricData.title
    if (metricData.subtitle) {
      title += `<br>${metricData.subtitle}`
    }

    const layout = {
      title,
      xaxis: {
        title: metricData.x.axis.title,
        zeroline: false
      },
      yaxis: {
        title: metricData.y.axis.title,
        zeroline: false
      }
    }

    if (metricData.type === 'bar' && metricData.y.data.length > 1) {
      layout.barmode = 'group'
      if (metricData.options && metricData.options.layout) {
        layout.barmode = metricData.options.layout
      }
    }

    if (metricData.x.axis.range) {
      layout.xaxis.range = metricData.x.axis.range
      if (!metricData.y.axis.range) {
        let yAxisRange = [0, 0]
        for (let i = 0; i < metricData.y.data.length; i += 1) {
          const yData = metricData.y.data[i].values
          const xData =
            metricData.x.data.length === metricData.y.data.length
              ? metricData.x.data[i].values
              : metricData.x.data[0].values
          for (const [x, y] of zip(xData, yData)) {
            if (
              x < metricData.x.axis.range[0] ||
              x > metricData.x.axis.range[1]
            ) {
              continue
            }
            if (y < yAxisRange[0]) {
              yAxisRange[0] = y
            }
            if (y > yAxisRange[1]) {
              yAxisRange[1] = y
            }
          }
        }
        layout.yaxis.range = yAxisRange
      }
    }

    if (metricData.y.axis.range) {
      layout.yaxis.range = metricData.y.axis.range
    }

    Plotly.newPlot(container, chartData, layout)
  }

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
                  return `<td title="${row[0]}">${
                    value === null ? '—' : value
                  }</td>`
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
                      `<td title="${tableData.data.columns[i]}">${
                        value === null ? '—' : value
                      }</td>`
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
    const isTransposed = tableElement && 'transposed' in tableElement.dataset
    summaryTable(summary, !isTransposed)
  }

  window.summaryDownload = summaryDownload
  function summaryDownload() {
    const tableElement = document.querySelector('#summary-table')
    const isTransposed = tableElement && 'transposed' in tableElement.dataset
    let dataArr
    if (isTransposed) {
      dataArr = zip(summary.data.columns, ...summary.data.rows)
    } else {
      dataArr = [summary.data.columns, ...summary.data.rows]
    }
    const csvData = csv.unparse(dataArr)
    const blob = new Blob([csvData], { type: 'text/plain;charset=utf-8' })
    saveAs(blob, 'alfred-summary-stats.csv')
  }

  function showExample(filename) {
    hideElement(resultContainer)
    hideElement(resultError)
    chartsContainer.innerHTML = ''
    summaryTab.innerHTML = ''
    showElement(resultInfo)
    resultLink.click()

    const url = `./examples/${filename}`
    axios
      .get(url, {
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

  function showError(message) {
    hideElement(resultContainer)
    hideElement(resultInfo)
    summaryTab.innerHTML = ''

    showElement(resultError)
    const errorMsg = document.querySelector('#error-message')
    if (errorMsg) errorMsg.innerHTML = message
  }

  function showElement(element) {
    if (element) element.classList.remove('d-none')
  }

  function hideElement(element) {
    if (element) element.classList.add('d-none')
  }
})
