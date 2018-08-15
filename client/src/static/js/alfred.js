import axios from 'axios'
import merge from 'deepmerge'
import pako from 'pako'

const API_URL = process.env.API_URL

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

let data, exampleData
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

function run() {
  const file = inputFile.files[0]

  // TODO error handling
  if (file === undefined) return

  hideElement(resultContainer)
  hideElement(resultError)
  showElement(resultInfo)

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
  const readGroups = {}
  data.samples.forEach(sample => {
    readGroups[sample.id] = sample.readGroups.map(rg => rg.id)
  })

  selectSample.innerHTML = samples.map(s => `<option>${s}</option>`).join('')

  selectReadGroup.innerHTML = readGroups[samples[0]]
    .map(rg => `<option>${rg}</option>`)
    .join('')

  vis(data, samples[0], readGroups[samples[0]][0])
}

window.handleSelectChange = handleSelectChange
function handleSelectChange() {
  const sample = selectSample.value
  const readGroup = selectReadGroup.value
  chartsContainer.innerHTML = ''
  vis(data, sample, readGroup)
}

const chartDispatch = {
  'line': lineChart,
  'table': table
}

function vis(data, sample, readGroup) {
  console.log(sample, readGroup)
  const dataRg = data.samples
    .find(s => s.id === sample)
    .readGroups.find(rg => rg.id === readGroup)

  for (const metric of dataRg.metrics) {
    console.log(metric.id)
  }
}

function lineChart(metricData) {
  const container = document.createElement('div')
  chartsContainer.appendChild(container)

  const xData = metricData.x.data[0]
  const chartData = []
  for (const yData of metricData.y.data) {
    const trace = {/*...*/}
    chartData.push(trace)
  }

  const layout = {
    title: 'TODO',
    xaxis: {
      title: 'TODO',
      zeroline: false
    },
    yaxis: {
      title: 'TODO',
      zeroline: false
    }
  }

  if (metricData.x.axis.range) {
    layout.xaxis.range = metricData.x.axis.range
  }

  if (metricData.y.axis.range) {
    layout.yaxis.range = metricData.y.axis.range
  }

  // Plotly.newPlot(container, chartData, layout)

  // const chartData = [{ x, y }]
  // const chartLayout = merge(
  //   {
  //     xaxis: {
  //       zeroline: false
  //     },
  //     yaxis: {
  //       zeroline: false
  //     }
  //   },
  //   layout
  // )
  // Plotly.newPlot(container, chartData, chartLayout)
}

function table() {
  console.log('TODO: summary table')
}

function renderBaseContentChart(data) {
  const container = document.createElement('div')
  chartsContainer.appendChild(container)
  const title = data.name
  const x = data.pos
  const chartData = []
  for (const base of 'ACGTN') {
    chartData.push({
      name: base,
      mode: 'scatter',
      x,
      y: data[base]
    })
  }
  const layout = {
    title,
    xaxis: {
      title: 'Position in read',
      zeroline: false
    },
    yaxis: {
      title: 'Base count',
      zeroline: false
    }
  }
  Plotly.newPlot(container, chartData, layout)
}

function lineChart(x, y, layout = {}) {
  const container = document.createElement('div')
  chartsContainer.appendChild(container)
  const chartData = [{ x, y }]
  const chartLayout = merge(
    {
      xaxis: {
        zeroline: false
      },
      yaxis: {
        zeroline: false
      }
    },
    layout
  )
  Plotly.newPlot(container, chartData, chartLayout)
}

function renderReadLengthChart(data) {
  lineChart(data.length, data.count, {
    title: data.name,
    xaxis: {
      title: 'Read length'
    },
    yaxis: {
      title: 'Count'
    }
  })
}

function renderBaseQualityChart(data) {
  lineChart(data.pos, data.qual, {
    title: data.name,
    xaxis: {
      title: 'Position in read'
    },
    yaxis: {
      title: 'Mean base quality'
    }
  })
}

function renderMappingQualityChart(data) {
  lineChart(data.pos, data.qual, {
    title: data.name,
    xaxis: {
      title: 'Mapping quality'
    },
    yaxis: {
      title: 'Count'
    }
  })
}

function renderCoverageChart(data) {
  lineChart(data.coverage, data.count, {
    title: data.name,
    xaxis: {
      title: 'Coverage',
      range: [0, 60]
    },
    yaxis: {
      title: 'Count'
    }
  })
}

function renderInsertSizeChart(data) {
  const container = document.createElement('div')
  chartsContainer.appendChild(container)
  const title = data.name
  const x = data.insertSize
  const chartData = []
  const orientationLabel = {
    fMinus: 'F-',
    fPlus: 'F+',
    rMinus: 'R-',
    rPlus: 'R+'
  }
  for (const orientation of ['fMinus', 'fPlus', 'rMinus', 'rPlus']) {
    chartData.push({
      name: orientationLabel[orientation],
      mode: 'scatter',
      x,
      y: data[orientation]
    })
  }
  const layout = {
    title,
    xaxis: {
      title: 'Insert size',
      range: [0, 1000],
      zeroline: false
    },
    yaxis: {
      title: 'Count',
      zeroline: false
    }
  }
  Plotly.newPlot(container, chartData, layout)
}

function renderOnTargetRateChart(data) {
  lineChart(data.targetExtension, data.fractionOnTarget, {
    title: data.name,
    xaxis: {
      title: 'Left/right extension of target region'
    },
    yaxis: {
      title: 'Fraction of reads on-target'
    }
  })
}

function renderTargetCoverageChart(data) {
  lineChart(data.coverageLevel, data.fractionAboveCoverage, {
    title: data.name,
    xaxis: {
      title: 'Coverage'
    },
    yaxis: {
      title: 'Fraction of targets above given coverage'
    }
  })
}

function showExample() {
  hideElement(resultContainer)
  chartsContainer.innerHTML = ''
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
