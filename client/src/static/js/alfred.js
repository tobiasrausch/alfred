import axios from 'axios'
import pako from 'pako'

$('#mainTab a').on('click', function(e) {
  e.preventDefault()
  $(this).tab('show')
})

const resultLink = document.getElementById('link-results')

const submitButton = document.getElementById('btn-submit')
submitButton.addEventListener('click', function() {
  resultLink.click()
  run()
})

let data, exampleData, readGroups
const exampleButton = document.getElementById('btn-example')
exampleButton.addEventListener('click', showExample)
const urlExample = document.getElementById('link-example')

const inputFile = document.getElementById('inputFile')
const linkPdf = document.getElementById('link-pdf')
const chartsContainer = document.getElementById('charts-container')
const resultContainer = document.getElementById('result-container')
const resultInfo = document.getElementById('result-info')
const resultError = document.getElementById('result-error')
const selectSample = document.getElementById('select-sample')
const selectReadGroup = document.getElementById('select-rg')
const summaryTab = document.getElementById('summary-tab')

function run() {
  const file = inputFile.files[0]

  // TODO error handling
  if (file === undefined) return

  hideElement(resultContainer)
  hideElement(resultError)
  showElement(resultInfo)
  summaryTab.innerHTML = ''

  const fileReader = new FileReader()
  const isGzip = /\.gz$/.test(file.name)
  if (isGzip) {
    fileReader.readAsArrayBuffer(file)
  } else {
    fileReader.readAsText(file)
  }
  fileReader.onload = event => {
    let content = event.target.result
    if (isGzip) {
      content = pako.ungzip(content, { to: 'string' })
    }
    data = JSON.parse(content)
    handleSuccess(data)
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

  table(
    {
      title: data.samples[0].summary.title,
      data: {
        columns: data.samples[0].summary.data.columns,
        rows: data.samples
          .map(s => s.summary.data.rows)
          .reduce((acc, cur) => acc.concat(cur))
      }
    },
    summaryTab
  )

  vis(data, samples[0], readGroups[samples[0]][0])
}

window.handleReadGroupSelectChange = handleReadGroupSelectChange
function handleReadGroupSelectChange() {
  const sample = selectSample.value
  const readGroup = selectReadGroup.value
  chartsContainer.innerHTML = ''
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
  vis(data, sample, readGroup)
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
  element.innerHTML = html
  parent.appendChild(element)
}

function showExample() {
  hideElement(resultContainer)
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

function showElement(element) {
  element.classList.remove('d-none')
}

function hideElement(element) {
  element.classList.add('d-none')
}