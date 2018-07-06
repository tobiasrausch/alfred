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

const exampleButton = document.getElementById('btn-example')
exampleButton.addEventListener('click', showExample)
const inputFile = document.getElementById('inputFile')

function run() {
  const formData = new FormData()
  formData.append('queryFile', inputFile.files[0])
  hideElement(resultContainer)
  hideElement(resultError)
  showElement(resultInfo)

  axios
    .post(`${API_URL}/upload`, formData)
    .then(res => {
      if (res.status === 200) {
        handleSuccess(`${API_URL}/` + res.data.data.url)
      }
    })
    .catch(err => {
      let errorMessage = err
      if (err.response) {
        errorMessage = err.response.data.errors
          .map(error => error.title)
          .join('; ')
      }
      hideElement(resultInfo)
      showElement(resultError)
      resultError.querySelector('#error-message').textContent = errorMessage
    })
}

async function handleSuccess(pdfUrl) {
  hideElement(resultInfo)
  hideElement(resultError)
  showElement(resultContainer)

  canvasContainer.innerHTML = ''

  linkPdf.href = pdfUrl
  const w = canvasContainer.clientWidth
  const pdf = await pdfjsLib.getDocument({ url: pdfUrl })

  for (let i = 1; i <= pdf.numPages; i += 1) {
    const page = await pdf.getPage(i)
    const canvas = document.createElement('canvas')
    canvas.width = w

    const scale = canvas.width / page.getViewport(1).width
    const viewport = page.getViewport(scale)
    canvas.height = page.getViewport(1).height * scale
    canvasContainer.appendChild(canvas)

    const context = canvas.getContext('2d')

    const renderContext = {
      canvasContext: context,
      viewport: viewport
    }

    page.render(renderContext)
  }
}

function showExample() {
  resultLink.click()
  handleSuccess(linkExample.href)
}

function showElement(element) {
  element.classList.remove('d-none')
}

function hideElement(element) {
  element.classList.add('d-none')
}
