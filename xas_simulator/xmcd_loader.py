"""
Load experimental data from beamlines

By Dan Porter
Diamond Light Source Ltd
2024
"""

import os
import re

import h5py
import numpy as np
from hdfmap import create_nexus_map

from nexpy.gui.datadialogs import GridParameters, NXDialog
from nexpy.gui.pyqt import QtWidgets
from nexpy.gui.utils import report_error, display_message, keep_data, nxload
from nexusformat.nexus import NeXusError, NXentry, NXdata, NXfield
from nexpy.gui.consoleapp import _tree

# from .xmcd_analysis import calculate_xmcd
from .beamline_data_loader import get_scan_paths
from .spectra_analysis import get_spectra, fit_bkg_then_norm_to_jump, fit_bkg_then_norm_to_peak, calc_subtraction


# DEFAULT_DATA_PATH = '/dls/i06-1/data/2024/'
DEFAULT_DATA_PATH = '/dls/science/groups/das/ExampleData/xmcd'
DEFAULT_PATH_SPEC = 'i06-1-%d.nxs'
EN_PATH = '/entry/instrument/fesData/pgmenergy'
SIGAL_PATH = '/entry/instrument/fesData/C1'
MONITOR_PATH = '/entry/instrument/fesData/C2'
POL_PATH = '/entry/instrument/id/polarisation'
PARAM_STRING = 'pol={polarisation}, T={T_sample:5.1f}K, B={field_sum:6.1f}T'


def show_dialog():
    try:
        dialog = DataLoader()
        dialog.show()
    except NeXusError as error:
        report_error("XAS data loader", error)


class DataLoader(NXDialog):

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

        # XMCD parameters
        self.xmcd = GridParameters()
        self.xmcd.add('first', 339089, 'First', spinbox=True)
        self.xmcd.add('last', 339093, 'Last', spinbox=True)
        self.xmcd.add('step', 1, 'Step', spinbox=True)
        # self.xmcd.add('energy_path', EN_PATH, 'Energy')
        self.xmcd.add('signal_path', SIGAL_PATH, 'Signal')
        # self.xmcd.add('monitor_path', MONITOR_PATH, 'Monitor')
        # self.xmcd.add('pol_path', POL_PATH, 'Polarisation')

        # action_buttons = self.action_buttons(('Get Ei', self.get_ei))
        self.set_layout(
            self.action_buttons(('Select Folder', self.load_directory)),
            self.experiment.grid(header=False, width=400),
            'stretch',
            self.action_buttons(('Load File', self.load_data)),
            self.xmcd.grid(header=False, width=300),
            self.action_buttons(('Calculate XMCD', self.calculate_xmcd))
        )
        self.set_title('XMCD Data Loader')

    def _msg(self, message):
        display_message('xmcd loader', message)

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
        self.xmcd['first'].value = last_scan - 10
        self.xmcd['last'].value = last_scan

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
        first = int(self.xmcd['first'].value)
        last = int(self.xmcd['last'].value)
        step = int(self.xmcd['step'].value)
        scan_range = range(first, last + step, step)
        scan_files = [os.path.join(path, spec % scanno) for scanno in scan_range]

        # energy_path = self.xmcd['energy_path'].value
        signal_path = self.xmcd['signal_path'].value
        # monitor_path = self.xmcd['monitor_path'].value
        # pol_path = self.xmcd['pol_path'].value
        #
        # av_energy, interp_pc, interp_nc, diff = calculate_xmcd(
        #     file_list=scan_files,
        #     energy_path=energy_path,
        #     signal_path=signal_path,
        #     monitor_path=monitor_path,
        #     pol_path=pol_path
        # )

        signal_name = os.path.basename(signal_path)  # default=None, or set to e.g. 'C1'
        energy_path, signal_path, monitor_path, pol_path = get_scan_paths(scan_files[0], signal_name)

        # load data, average and normalise spectra
        av_energy, signals, raw = get_spectra(
            *scan_files,
            energy_path=energy_path,
            signal_path=signal_path,
            pol_path=pol_path,
            monitor_path=monitor_path
        )

        # normalise spectra and remove background (normalise high energy region to 1)
        norm_signals = {
            # label: fit_bkg_then_norm_to_peak(
            label: fit_bkg_then_norm_to_jump(
                energy=av_energy,
                signal=signal,
                ev_from_start=6,  # defines background region below edge, in eV from start of scan
                ev_from_end=10  # defines high-energy region above edge, in eV from end of scan
            )[0] for label, signal in signals.items()
        }

        # determine polarisations and calculate subtraction (XMCD or XMLD)
        result, subtraction_type = calc_subtraction(norm_signals)



        # Load files for tree
        pc = NXentry(name='pc')
        nc = NXentry(name='nc')
        # lv = NXentry(name='lv')
        # lh = NXentry(name='lh')
        regex = re.compile(spec.replace('%d', '(\d+)'))
        for file in scan_files:
            scan_number = regex.findall(os.path.basename(file))[0]
            nx = nxload(file)
            pol = nx[pol_path]
            if pol == 'pc':
                pc[scan_number] = nx['entry']
            elif pol == 'nc':
                nc[scan_number] = nx['entry']
            # elif pol == 'lv':
            #     lv[scan_number] = nx['entry']
            # elif pol == 'lh':
            #     lh[scan_number] = nx['entry']
            else:
                self._msg(f"file {scan_number} has an unknown polarisation: {pol}")

        # for plotting
        en = NXfield(av_energy, name='Energy_eV')
        xas1 = NXdata(NXfield(signals['pc'], name='XAS_pc'), (en,), name='XAS_pc')
        xas2 = NXdata(NXfield(signals['nc'], name='XAS_nc'), (en,), name='XAS_nc')
        xmcd = NXdata(NXfield(result, name=subtraction_type), (en,), name=subtraction_type)
        entry = NXentry(pc, nc, xas1, xas2, xmcd, name='processed')
        xmcd.set_default()

        # entry.set_default()
        xas1.plot('b-', label='pc')
        xas2.oplot('g-', label='nc')
        xmcd.oplot('r-', label='XMCD')

        keep_data(entry)


