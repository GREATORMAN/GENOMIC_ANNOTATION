from flask import Flask, render_template, request, redirect, url_for, flash
import subprocess
import os
from Bio import SeqIO

app = Flask(__name__)
app.secret_key = "supersecretkey"

UPLOAD_FOLDER = 'uploads'
os.makedirs(UPLOAD_FOLDER, exist_ok=True)

bioseq_path = r"C:\Strawberry\perl\site\bin\bioseq.bat"


def format_bioseq_output(raw_output):
    lines = raw_output.strip().split('\n')
    rows = [line.split() for line in lines if line.strip()]
    html = "<table class='table table-striped table-bordered'>"
    for row in rows:
        html += "<tr>" + "".join(f"<td>{col}</td>" for col in row) + "</tr>"
    html += "</table>"
    return html


def format_features_output(features):
    html = ""
    for feat in features:
        # Split the feature string by <br> to separate key:value pairs
        parts = feat.split('<br>')
        html += '<div class="gene-block">'
        for part in parts:
            if ':' in part:
                key, val = part.split(':', 1)
                html += f'<p><strong>{key.strip()}:</strong> {val.strip()}</p>'
            else:
                # In case line does not contain ':', just add as is
                html += f'<p>{part.strip()}</p>'
        html += '</div>'
    return html


@app.route('/')
def home():
    return render_template('index.html')


@app.route('/analyze', methods=['POST'])
def analyze():
    file = request.files.get('fasta_file')
    operation = request.form.get('operation')
    subseq_range = request.form.get('subseq_range')
    reloop_pos = request.form.get('reloop_pos')

    if not file or file.filename == '':
        flash('No file selected!', 'danger')
        return redirect(url_for('home'))

    if not operation:
        flash('No operation selected!', 'danger')
        return redirect(url_for('home'))

    filepath = os.path.join(UPLOAD_FOLDER, file.filename)
    file.save(filepath)

    if operation == 'feat2fas':
        try:
            record = SeqIO.read(filepath, "genbank")
            features = []
            for feat in record.features:
                if feat.type == "gene":
                    loc = str(feat.location)
                    strand = feat.strand
                    qualifiers = feat.qualifiers
                    gene_name = qualifiers.get('gene', ['Unnamed'])[0]
                    features.append(
                        f"Gene: {gene_name}<br>"
                        f"Location: {loc}<br>"
                        f"Strand: {strand}<br>"
                        f"Qualifiers: {qualifiers}<br>"
                    )
            result = format_features_output(features) if features else "No gene features found in the GenBank file."
        except Exception as e:
            result = f"Error reading GenBank file: {e}"
        return render_template('result.html', output=result, filename=file.filename)

    command = ['cmd', '/c', bioseq_path]

    if operation == 'gc':
        command += ['-c', filepath]
    elif operation == 'length':
        command += ['-l', filepath]
    elif operation == 'count':
        command += ['-n', filepath]
    elif operation == 'revcom':
        command += ['-r', filepath]
    elif operation == 'translate1':
        command += ['-t', '1', filepath]
    elif operation == 'translate3':
        command += ['-t', '3', filepath]
    elif operation == 'translate6':
        command += ['-t', '6', filepath]
    elif operation == 'mol_wt':
        command += ['--mol-wt', filepath]
    elif operation == 'iep':
        command += ['--iep', filepath]
    elif operation == 'lead_gaps':
        command += ['--lead-gaps', filepath]
    elif operation == 'gaps_dna':
        command += ['--num-gaps-dna', filepath]
    elif operation == 'gaps_aa':
        command += ['--num-gaps-aa', filepath]
    elif operation == 'subseq':
        if not subseq_range or ',' not in subseq_range:
            flash("Please provide start,end for subsequence.", 'danger')
            return redirect(url_for('home'))
        command += ['-s', subseq_range, filepath]
    elif operation == 'remove_gaps':
        command += ['--no-gaps', filepath]
    elif operation == 'remove_stop':
        command += ['--remove-stop', filepath]
    elif operation == 'reloop':
        if not reloop_pos:
            flash("Please provide a position for relooping.", 'danger')
            return redirect(url_for('home'))
        command += ['--reloop', reloop_pos, filepath]
    elif operation == 'linearize':
        command += ['--linearize', filepath]
    elif operation == 'longest_orf':
        command += ['--longest-orf', filepath]
    else:
        flash('Unsupported operation selected!', 'danger')
        return redirect(url_for('home'))

    try:
        raw_result = subprocess.check_output(command, stderr=subprocess.STDOUT, shell=False, timeout=20).decode('utf-8')
        result = format_bioseq_output(raw_result)
    except subprocess.CalledProcessError as e:
        result = f"Error running bioseq:\n{e.output.decode('utf-8')}"
    except Exception as e:
        result = f"Unexpected error: {str(e)}"

    return render_template('result.html', output=result, filename=file.filename)


@app.route('/genome-annotate', methods=['GET', 'POST'])
def genome_annotate():
    if request.method == 'POST':
        file = request.files.get('fasta_file')
        if not file or file.filename == '':
            flash('No file selected!', 'danger')
            return redirect(url_for('genome_annotate'))

        filepath = os.path.join(UPLOAD_FOLDER, file.filename)
        file.save(filepath)

        try:
            record = SeqIO.read(filepath, "genbank")
            features = []
            for feat in record.features:
                if feat.type in ("gene", "CDS"):
                    loc = str(feat.location)
                    strand = feat.strand
                    qualifiers = feat.qualifiers
                    gene_name = qualifiers.get('gene', ['Unnamed'])[0]
                    product = qualifiers.get('product', [''])[0]
                    protein_id = qualifiers.get('protein_id', [''])[0]
                    note = qualifiers.get('note', [''])[0]
                    features.append(
                        f"Feature type: {feat.type}<br>"
                        f"Gene: {gene_name}<br>"
                        f"Location: {loc}<br>"
                        f"Strand: {strand}<br>"
                        f"Product: {product}<br>"
                        f"Protein ID: {protein_id}<br>"
                        f"Note: {note}<br>"
                        f"Qualifiers: {qualifiers}<br>"
                    )
            result = format_features_output(features) if features else "No gene or CDS features found in the GenBank file."
        except Exception as e:
            result = f"Error reading GenBank file: {e}"

        return render_template('genome_annotate_result.html', output=result, filename=file.filename)

    return render_template('genome_annotate.html')


if __name__ == '__main__':
    app.run(debug=True)
