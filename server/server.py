import os
import uuid
import re
import subprocess
from subprocess import call
from flask import Flask, send_file, flash, render_template, request, redirect, url_for
from werkzeug.utils import secure_filename


ALFREDWS = os.path.dirname(os.path.abspath(__file__))

app = Flask(__name__)
app.config['ALFRED'] = os.path.join(ALFREDWS, "..")
app.config['UPLOAD_FOLDER'] = os.path.join(app.config['ALFRED'], "data")
app.config['MAX_CONTENT_LENGTH'] = 8 * 1024 * 1024   #maximum of 8MB
app.secret_key = 'soadfdafvmv'

def allowed_file(filename):
   return '.' in filename and filename.rsplit('.', 1)[1].lower() in set(['pdf', 'gz'])

uuid_re = re.compile(r'^[0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12}$')
def is_valid_uuid(s):
   return uuid_re.match(s) is not None

@app.route('/download/<uuid>')
def download(uuid):
   if is_valid_uuid(uuid):
      filename = "alfred_" + uuid + ".pdf"
      if allowed_file(filename):
         sf = os.path.join(app.config['UPLOAD_FOLDER'], uuid[0:2])
         if os.path.exists(sf):
            if os.path.isfile(os.path.join(sf, filename)):
               return send_file(os.path.join(sf, filename), attachment_filename=filename)
   return "File does not exist!"

@app.route('/upload', methods = ['GET', 'POST'])
def upload_file():
   if request.method == 'POST':
      uuidstr = str(uuid.uuid4())

      # Get subfolder
      sf = os.path.join(app.config['UPLOAD_FOLDER'], uuidstr[0:2])
      if not os.path.exists(sf):
         os.makedirs(sf)

      # Statistics file
      if 'stats' not in request.files:
         error = "Alfred statistic file missing!"
         return render_template('upload.html', error = error)
      fexp = request.files['stats']
      if fexp.filename == '':
         error = "Alfred statistic file missing!"
         return render_template('upload.html', error = error)
      if not allowed_file(fexp.filename):
         error = "Alfred statistic file has incorrect file type!"
         return render_template('upload.html', error = error)
      fexpname = os.path.join(sf, "alfred_" + uuidstr + "_" + secure_filename(fexp.filename))
      fexp.save(fexpname)

      # Run Rscript
      outfile = os.path.join(sf, "alfred_" + uuidstr + ".pdf")
      logfile = os.path.join(sf, "alfred_" + uuidstr + ".log")
      errfile = os.path.join(sf, "alfred_" + uuidstr + ".err")
      with open(logfile, "w") as log:
         with open(errfile, "w") as err:
            blexe = os.path.join(app.config['ALFRED'], "R/stats.R")
            return_code = call(['Rscript', blexe, fexpname, outfile], stdout=log, stderr=err)
      if return_code != 0:
         error = "Error in running Alfred!"
         return render_template('upload.html', error = error)

      # Send download pdf
      return redirect("/alfred/download/" + uuidstr, code=302)
   return render_template('upload.html')

@app.route("/")
def submit():
    return render_template("upload.html")

if __name__ == '__main__':
   app.run(host = '0.0.0.0', port = 3300, debug = True, threaded=True)
