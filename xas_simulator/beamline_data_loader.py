"""
Load data from beamlines
"""

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
    'energy': '(fastEnergy|pgm_energy|energy)',
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


def get_scan_paths(filename):
    """
    Return hdf dataset paths
    """
    nxmap = hdfmap.create_nexus_map(filename)

    energy_path = nxmap.scannables['pgmenergy']
    signal_path = nxmap.scannables['C3']
    monitor_path = nxmap.scannables['C4']
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

