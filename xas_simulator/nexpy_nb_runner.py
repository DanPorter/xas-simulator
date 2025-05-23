"""
Nexpy widget for notebook runner
"""

import os
import re

from nexpy.gui.datadialogs import GridParameters, NXDialog
from nexpy.gui.pyqt import QtWidgets
from nexpy.gui.utils import report_error, display_message
from nexusformat.nexus import NeXusError

from xas_simulator.notebook_runner import xas_notebook, launch_server, launch_browser

# DEFAULT_DATA_PATH = '/dls/i06-1/data/2024/'
# DEFAULT_DATA_PATH = '/scratch/grp66007/data/i06/mm36794-1'
DEFAULT_DATA_PATH = '/dls/science/groups/das/ExampleData/xmcd'
DEFAULT_PATH_SPEC = 'i06-1-%d.nxs'
# DEFAULT_OUTPUT = '/scratch/grp66007/data/i06/mm36794-1/processing'
DEFAULT_OUTPUT = '/dls/science/groups/das/ExampleData/xmcd'
DEFAULT_QUANTY_PATH = 'Quanty'


def show_dialog():
    try:
        dialog = NBRunner()
        dialog.show()
    except NeXusError as error:
        report_error("NBRunner", error)


class NBRunner(NXDialog):

    def __init__(self, parent=None):

        super().__init__(parent)

        # File parameters
        self.experiment = GridParameters()
        self.experiment.add('Data path', DEFAULT_DATA_PATH, 'Experiment data path')
        self.experiment.add('path spec', DEFAULT_PATH_SPEC, 'Data file specifier')
        self.experiment.add('output', DEFAULT_OUTPUT, 'Output path')

        # Data file parameters
        self.files = GridParameters()
        self.files.add('first', 336801, 'First', spinbox=True)
        self.files.add('last', 336804, 'Last (inclusive)', spinbox=True)
        self.files.add('step', 1, 'Step', spinbox=True)
        self.files.add('detector_name', 'default', 'Detector Name')

        self.parameters = GridParameters()
        self.parameters.add('qtypath', DEFAULT_QUANTY_PATH, 'Quanty path')
        self.parameters.add('ion', 'Ni', 'Ion Name')
        self.parameters.add('charge', 2, 'Charge')
        self.parameters.add('sym', 'Oh', 'Symmetry')
        self.parameters.add('beta', 0.8, 'Beta')
        self.parameters.add('10Dq', 1.0, '10Dq')
        self.parameters.add('Bx', 0.0, 'B_x')
        self.parameters.add('By', 0.0, 'B_y')
        self.parameters.add('Bz', 0.0, 'B_z')
        self.parameters.add('Hx', 0.0, 'H_x')
        self.parameters.add('Hy', 0.0, 'H_y')
        self.parameters.add('Hz', 0.1, 'H_z')
        self.parameters.add('T', 1.0, 'T (K)')

        # action_buttons = self.action_buttons(('Get Ei', self.get_ei))
        self.set_layout(
            self.action_buttons(('Select Folder', self.load_directory)),
            self.experiment.grid(header=False, width=400),
            'stretch',
            self.files.grid(header=False, width=300),
            'stretch',
            self.parameters.grid(header=False, width=300),
            self.action_buttons(('Run Notebook', self.run_notebook)),
            self.action_buttons(('View Notebook', self.launch_notebook)),
            self.action_buttons(('Launch Jupyter', self.launch_jupyter))
        )
        self.set_title('XMCD Data Loader')

    def _msg(self, message):
        display_message('NBRunner', message)

    def load_directory(self):
        directory = self.experiment['Data path'].value
        directory = QtWidgets.QFileDialog.getExistingDirectory(
            self, 'Choose Directory', directory)
        if directory is None or not os.path.exists(directory):
            return
        self.experiment['Data path'].value = directory
        self.experiment['output'].value = directory + '/processing'
        last_scan = self.get_last_scan_number()
        self.files['first'].value = last_scan - 10
        self.files['last'].value = last_scan

    def get_last_scan_number(self) -> int:
        path = self.experiment['Data path'].value
        spec = self.experiment['path spec'].value
        files = os.listdir(path)
        file_str = '\n'.join(files)
        # convert path spec to regex
        regex = spec.replace('%d', '(\d+)')
        scan_numbers = re.findall(regex, file_str)
        if scan_numbers:
            scan_numbers.sort()
            return int(scan_numbers[-1])
        return 0

    def run_notebook(self):
        path = self.experiment['Data path'].value
        spec = self.experiment['path spec'].value
        output_path = self.experiment['output'].value
        first = int(self.files['first'].value)
        last = int(self.files['last'].value)
        step = int(self.files['step'].value)
        detector_name = str(self.files['detector_name'].value)
        # filenames = [
        #     path for scanno in range(first, last + step, step)
        #     if h5py.is_hdf5(path := os.path.join(path, spec % scanno))
        # ]
        filenames = [
            os.path.join(path, spec % scanno) for scanno in range(first, last + step, step)
        ]
        output_name = f'xas_notebook_{first}_{last}.ipynb'
        output_file = os.path.join(output_path, output_name)

        if os.path.isfile(output_file):
            response = self.confirm_action(f"File: '{output_file}' already exists and will be overwritten, continue?")
            if not response:
                return

        qtypath = self.parameters['qtypath'].value
        ion = self.parameters['ion'].value
        ch = int(self.parameters['charge'].value)
        sym = self.parameters['sym'].value
        beta = self.parameters['beta'].value
        Dq = self.parameters['10Dq'].value
        mag_field = [self.parameters['Bx'].value,
                     self.parameters['By'].value,
                     self.parameters['Bz'].value]
        exchange_field = [self.parameters['Hx'].value,
                          self.parameters['Hy'].value,
                          self.parameters['Hz'].value]
        temperature = self.parameters['T'].value
        simulation = {
            'ion': ion,
            'ch': ch,
            'sym': sym,
            'beta': beta,
            'dq': Dq,
            'mag_field': mag_field,
            'exchange_field': exchange_field,
            'temperature': temperature,
            'quanty_path': qtypath,
            'detector_name': detector_name,
        }
        xas_notebook(
            files=filenames,
            simulation=simulation,
            output_file=output_file
        )
        self.display_message(f'Calculation finished, created\n{output_file}')

    def launch_notebook(self):
        """
        Launch jupyter notebook server
        """
        path = self.experiment['Data path'].value
        spec = self.experiment['path spec'].value
        output_path = self.experiment['output'].value
        first = int(self.files['first'].value)
        last = int(self.files['last'].value)
        # output_name = f'xas_notebook_{first}_{last}.ipynb'
        output_name = f'xas_notebook_{first}_{last}.html'
        output_file = os.path.join(output_path, output_name)

        if os.path.isfile(output_file):
            # launch_server(output_file)
            launch_browser(output_file)
        else:
            self.display_message(f'Run the calculation first!')

    def launch_jupyter(self):
        """
        Launch jupyter notebook server in the processing path
        """
        output_path = self.experiment['output'].value
        launch_server(output_path)