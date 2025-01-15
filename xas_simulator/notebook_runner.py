"""
Jupyter notebook runner

 - copy specified notebook to new desination
 - replace the first cell of the notebook
 - execute the notebook
 - save a html file of the executed notebook
"""

import os
import subprocess
import nbformat
from nbconvert.preprocessors import ExecutePreprocessor
from nbconvert import HTMLExporter


def process_notebook(template_notebook: str, output_notebook: str, first_cell_code: str):
    """
    Copy a notebook, replace the first cell with new code, then run the notebook

    :param template_notebook: str, filename of the jupyter notebook (.ipynb) to copy
    :param output_notebook: str, filename of the jupyter notebook (.ipynb) to save
    :param first_cell_code: str, python code to replace in the 1st cell of the template notebook
    """

    # output files
    name, ext = os.path.splitext(output_notebook)
    output_notebook = name + '.ipynb'
    output_html = name + '.html'

    # Open old notebook
    print('Copying old notebook: %s' % template_notebook)
    nb = nbformat.read(template_notebook, as_version=4)

    # Replace first cell
    print('Replacing code')
    nb.cells[0]['source'] = first_cell_code

    # Create notebook exporter
    print('Processing...')
    processor = ExecutePreprocessor()
    html_exporter = HTMLExporter()
    processor.preprocess(nb, resources={'metadata': {}})
    # processor.preprocess(nb, resources={'metadata': {'path': 'notebooks/'}})  # specify working directory

    with open(output_notebook, 'w', encoding='utf-8') as f:
        nbformat.write(nb, f)

    # save html
    (body, resources) = html_exporter.from_notebook_node(nb)
    with open(output_html, 'w') as f:
        f.write(body)
    print('Completed. HTML written to %s' % output_html)


def launch_server(notebook_file = ''):
    """
    Launch jupyter notebook server for a single file
    """
    shell_cmd = f"gnome-terminal -- bash -c \"jupyter notebook {notebook_file}; exec bash\""
    subprocess.Popen(shell_cmd, shell=True)

def xas_notebook(files: list, simulation: dict, output_file: str):
    """
    process an xas notebook
    """
    cell_str = (
            "data_files = [\n    " +
            '\n'.join(f"    '{file}'," for file in files) +
            '\n]\n' +
            f"sim_data = {str(simulation)}"
    )
    notebook_path = os.path.join(os.path.dirname(__file__), 'notebooks')
    process_notebook(os.path.join(notebook_path, 'xas_notebook.ipynb'), output_file, cell_str)
