"""
Load experimental data from beamlines

By Dan Porter
Diamond Light Source Ltd
2024
"""

import os
import re

import h5py
import hdfmap
import numpy as np
from hdfmap import create_nexus_map

from nexpy.gui.datadialogs import GridParameters, NXDialog
from nexpy.gui.pyqt import QtWidgets
from nexpy.gui.utils import report_error, display_message, keep_data, nxload
from nexusformat.nexus import NeXusError, NXentry, NXdata, NXfield
from nexpy.gui.consoleapp import _tree

from .xmcd_analysis import calculate_xmcd


DEFAULT_DATA_PATH = '/dls/i06/data/2024/'
DEFAULT_PATH_SPEC = 'i06-%d.nxs'
PARAM_STRING = 'pol={polarisation}, {incident_energy:.2f} eV'


def show_dialog():
    try:
        dialog = PeemLoader()
        dialog.show()
    except NeXusError as error:
        report_error("PEEM data loader", error)


class PeemLoader(NXDialog):

    def __init__(self, parent=None):

        super().__init__(parent)

        self.nexus_map = None  # hdfmap.NexusMap

        # File parameters
        self.experiment = GridParameters()
        self.experiment.add('Data path', DEFAULT_DATA_PATH, 'Experiment data path')
        self.experiment.add('path spec', DEFAULT_PATH_SPEC, 'Data file specifier')
        self.experiment.add('scan number', 0, 'Scan number',
                            spinbox=True, slot=self.get_cmd)
        self.experiment.add('cmd', '', 'Scan command')
        self.experiment.add('params', '', 'Parameters')


        # action_buttons = self.action_buttons(('Get Ei', self.get_ei))
        self.set_layout(
            self.action_buttons(('Select Folder', self.load_directory)),
            self.experiment.grid(header=False, width=400),
            'stretch',
            self.action_buttons(('Load File', self.load_data)),
            self.action_buttons(('Calculate XMCD', self.calculate_xmcd))
        )
        self.set_title('PEEM Data Loader')

    def _msg(self, message):
        display_message('PEEM loader', message)

    def get_cmd(self):
        path = self.experiment['Data path'].value
        spec = self.experiment['path spec'].value
        scanno = self.experiment['scan number'].value
        scan_file = os.path.join(path, spec % scanno)
        cmd = ''
        params = ''
        if h5py.is_hdf5(scan_file):
            if self.nexus_map is None:
                self.nexus_map = create_nexus_map(scan_file)
            with h5py.File(scan_file, 'r') as hdf:
                cmd = self.nexus_map.get_data(hdf, 'scan_command')
                params = self.nexus_map.format_hdf(hdf, PARAM_STRING)
        self.experiment['cmd'].value = cmd
        self.experiment['params'].value = params

    def load_data(self):
        path = self.experiment['Data path'].value
        spec = self.experiment['path spec'].value
        scanno = self.experiment['scan number'].value
        scan_file = os.path.join(path, spec % scanno)
        if h5py.is_hdf5(scan_file):
            nx = _tree.load(scan_file)
        else:
            self._msg(f"{scan_file} is not a NeXus file!")

    def load_directory(self):
        directory = self.experiment['Data path'].value
        directory = QtWidgets.QFileDialog.getExistingDirectory(
            self, 'Choose Directory', directory)
        if directory is None or not os.path.exists(directory):
            return
        self.experiment['Data path'].value = directory
        last_scan = self.get_last_scan_number()
        self.experiment['scan number'].value = last_scan

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

    def calculate_xmcd(self):
        path = self.experiment['Data path'].value
        spec = self.experiment['path spec'].value
        scanno = self.experiment['scan number'].value
        scan_file = os.path.join(path, spec % scanno)

        if not h5py.is_hdf5(scan_file):
            self._msg(f"'{scan_file}' is not a file")

        scan = hdfmap.NexusLoader(scan_file)
        MEDIPIX_PATH = '/entry/medipix'
        if not MEDIPIX_PATH in scan.map.groups:
            self._msg(f"'{scan_file}' does not contain medipix data")

        data = scan(MEDIPIX_PATH + '/data')
        pol = scan(MEDIPIX_PATH + '/pol').astype(str)
        axes = scan.map.get_attr(MEDIPIX_PATH, 'axes').astype(str).tolist()
        # Find polarisation images
        nc = data[pol == 'nc']  # (len(ds) * len(energy), :, :)
        pc = data[pol == 'pc']
        xmcd = nc - pc
        # Subtract background (lower) energy
        if 'energy' in axes:
            energy = scan(MEDIPIX_PATH + '/energy')
            en_pol = energy[pol == 'pc']  # len(ds) * len(energy)
            xmcd_image = (np.sum(xmcd[en_pol > np.mean(energy), :, :], axis=0) -
                          np.sum(xmcd[en_pol < np.mean(energy), :, :], axis=0))
        else:
            xmcd_image = np.sum(xmcd, axis=0)

        # Load files for tree
        regex = re.compile(spec.replace('%d', '(\d+)'))
        scan_number = regex.findall(os.path.basename(scan_file))[0]
        # nx_data = nxload(scan_file)['entry']

        # for plotting
        # xas1 = NXdata(NXfield(nc, name='images_nc'), name='images_nc')
        # xas2 = NXdata(NXfield(pc, name='images_pc'), name='images_pc')
        xmcd = NXdata(NXfield(xmcd_image, name='XMCD'), name='XMCD')
        # entry = NXentry(nx_data, xas1, xas2, xmcd, name=f'XMCD_{scan_number}')
        entry = NXentry(xmcd, name=f'XMCD_{scan_number}')
        xmcd.set_default()

        # entry.set_default()
        # xmcd.implot()

        keep_data(entry)

