"""
Load data from beamlines
"""

import os
import numpy as np
import matplotlib.pyplot as plt

from IPython.display import display, Markdown

import hdfmap


# Define metadata using hdfmap namespace commands
# see hdfmap.eval_functions.eval_hdf
VALUES = {
    'cmd': '(cmd|scan_command)',
    'date': 'start_time.strftime("%a-%d/%b/%Y %H:%M:%S")',
    'pol': 'polarisation?("lh")',
    'iddgap': 'iddgap',
    'rowphase': 'idutrp if iddgap == 100 else iddtrp',
    'endstation': 'instrument_name',
    'temp': '(T_sample|lakeshore336_sample?(300))',
    'rot': '(scmth|xabs_theta|ddiff_theta?(0))',
    'field': 'field_x?(0), field_y?(0), field_z?(0)',
    'energy': '(fastEnergy|pgm_energy|pgmenergy|energy)',
    'monitor': '(C2|ca62sr)',
    'tey': '(C1|ca61sr) / (C2|ca62sr)',
    'tfy': '(C3|ca63sr) / (C2|ca62sr)',
}

LABELS = {
    'scan': (
        '{filename} pol={polarisation?("lh")}, ' +
        'T={(T_sample|lakeshore336_sample?(300)):5.1f}K, ' +
        'B={np.sqrt(np.sum(np.square([field_x?(0), field_y?(0), field_z?(0)]))):6.1f}T, {scan_command}'
    ),
    'title': (
        '{filename} pol={polarisation?("lh")}, ' +
        'T={(T_sample|lakeshore336_sample?(300)):5.1f}K, ' +
        'B={np.sqrt(np.sum(np.square([field_x?(0), field_y?(0), field_z?(0)]))):6.1f}T\n{scan_command}'
    ),
    'label': '{filename}: {scan_command}',
    'str_temp': 'T = {(T_sample|lakeshore336_sample?(300)):.2f} K',
    'str_field': 'B = {np.sqrt(np.sum(np.square([field_x?(0), field_y?(0), field_z?(0)]))):.2f} T',
    'str_energy': 'E = {np.mean((fastEnergy|pgm_energy|energy?(0))):.2f} eV',
}

MARKDOWN = """
### {filename}
{start_time.strftime("%a-%d/%b/%Y %H:%M:%S")}

 - cmd = *{scan_command}*
 - pol = {polarisation?("lh")}
 - T = {(T_sample|lakeshore336_sample?(300)):.2f} K
 - B = {np.sqrt(np.sum(np.square([field_x?(0), field_y?(0), field_z?(0)]))):.2f} T
 - E = {np.mean((fastEnergy|pgm_energy|pgmenergy|energy?(0))):.2f} eV
"""


def get_scan_paths(filename, signal_name=None):
    """
    Return hdf dataset paths
    """
    nxmap = hdfmap.create_nexus_map(filename)

    if 'pgmenergy' in nxmap.scannables: # i06-1
        signal_name = signal_name if signal_name in nxmap.scannables else 'C1'
        energy_path = nxmap.scannables['pgmenergy']
        signal_path = nxmap.scannables[signal_name]  # TEY
        monitor_path = nxmap.scannables['C2']
    elif 'pgm_energy' in nxmap.scannables: # i10-1
        signal_name = signal_name if signal_name in nxmap.scannables else 'mcse17'
        energy_path = nxmap.scannables['pgm_energy']
        signal_path = nxmap.scannables[signal_name]  # TEY
        monitor_path = nxmap.scannables['mcse16']
    else:
        raise Exception('currently only configured for i06-1 and i10-1')
    pol_path = nxmap['polarisation']
    return energy_path, signal_path, monitor_path, pol_path


def display_scans(*files):
    """
    Return string of file details
    """
    return '\n'.join(hdfmap.hdf_format(list(files), LABELS['scan']))


def get_data(filename: str, nxmap: hdfmap.NexusMap = None) -> dict:
    """
    Return data from file
    """
    if nxmap is None:
        nxmap = hdfmap.create_nexus_map(filename)

    with nxmap.load_hdf() as nxs:
        values = {name: nxmap.eval(nxs, expr) for name, expr in VALUES.items()}
        labels = {name: nxmap.format_hdf(nxs, expr) for name, expr in LABELS.items()}
    return {**values, **labels}


class Experiment:
    def __init__(self, data_directory: str, filespec: str = '%d.nxs'):
        if not os.path.isdir(data_directory):
            raise FileNotFoundError(f"{data_directory} does not exist")
        self.data_directory = data_directory
        self.filespec = filespec
        self.markdown_format = MARKDOWN

    def all_files(self):
        return [os.path.basename(path) for path in hdfmap.list_files(self.data_directory)]

    def get_file(self, filename: str | None = None, scan_number: int | None = None) -> str:
        if scan_number:
            filename = self.filespec % scan_number
        return os.path.join(self.data_directory, filename)

    def read_scan(self, filename: str):
        filepath = self.get_file(filename)

        nxmap = hdfmap.create_nexus_map(filepath)

        axes_paths, signal_paths = nxmap.nexus_defaults()
        axes_name = nxmap.datasets[axes_paths[0]].name
        signal_name = nxmap.datasets[signal_paths[0]].name

        with nxmap.load_hdf() as hdf:
            axes_data = nxmap.eval(hdf, axes_name)
            signal_data = nxmap.eval(hdf, signal_name)
            scan_data = nxmap.get_scannables(hdf)
            rsp = nxmap.format_hdf(hdf, self.markdown_format)
            energy = nxmap.eval(hdf, '(pgmenergy|pgm_energy)')
            xas = nxmap.eval(hdf, '(C1|mcse17) / (C2|mcse16)')

        if not issubclass(type(axes_data), np.ndarray) or not np.issubdtype(axes_data.dtype, np.number) or len(
                axes_data) != nxmap.scannables_length():
            axes_name = '![fail]' + axes_name
            axes_data = np.arange(nxmap.scannables_length())
        if not issubclass(type(signal_data), np.ndarray) or not np.issubdtype(signal_data.dtype, np.number) or len(
                signal_data) != nxmap.scannables_length():
            signal_name = '![fail]' + signal_name
            signal_data = np.ones(nxmap.scannables_length())

        fig, [ax1, ax2] = plt.subplots(1, 2, figsize=[10, 4], dpi=60)
        ax1.plot(axes_data, signal_data, '-o')
        ax1.set_xlabel(axes_name)
        ax1.set_ylabel(signal_name)

        ax2.plot(energy, xas, '-')
        ax2.set_xlabel('Energy [eV]')
        ax2.set_ylabel('TEY / Monitor')

        plt.show()
        display(Markdown(rsp))