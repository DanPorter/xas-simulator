"""
XMCD Analysis
"""

import numpy as np
from scipy.optimize import minimize
import h5py


def average_energy_scans(*args: tuple[np.ndarray]):
    """Return the minimum range covered by all input arguments"""
    min_energy = np.max([np.min(en) for en in args])
    max_energy = np.min([np.max(en) for en in args])
    min_step = np.min([np.min(np.abs(np.diff(en))) for en in args])
    return np.arange(min_energy, max_energy + min_step, min_step)


def combine_energy_scans(energy, *args: tuple[np.ndarray, np.ndarray]):
    """Average energy scans, interpolating at given energy"""
    data = np.zeros([len(args), len(energy)])
    for n, (en, dat) in enumerate(args):
        data[n, :] = np.interp(energy, en, dat)
    return data.mean(axis=0)


def calculate_xmcd(file_list: list[str], energy_path: str, signal_path: str, pol_path: str, monitor_path: str = None):
    """Combine pc + nc scans separately, then subtract them"""
    en_list = []
    pc = []
    nc = []
    for file in file_list:
        with h5py.File(file, 'r', swmr=True) as hdf:
            energy = hdf[energy_path][()]
            signal = hdf[signal_path][()]
            if monitor_path:
                monitor = hdf[monitor_path][()]
                signal = signal / monitor
            pol = hdf[pol_path].asstr()[()]

        en_list.append(energy)
        if pol == 'pc':
            pc.append((energy, signal))
        else:
            nc.append((energy, signal))
    av_energy = average_energy_scans(*en_list)
    interp_pc = combine_energy_scans(av_energy, *pc)
    interp_nc = combine_energy_scans(av_energy, *nc)
    diff = interp_nc - interp_pc
    return av_energy, interp_pc, interp_nc, diff


def spectra_difference(energy_1, signal_1, energy_2, signal_2):
    """Compare 2 spectra, return difference/steps"""
    if min(energy_2) > max(energy_1) or max(energy_2) < min(energy_1):
        raise Exception(f"Energies dont match: ({min(energy_1)}, {max(energy_1)}), ({min(energy_2)}, {max(energy_2)})")
    energy = average_energy_scans(energy_1, energy_2)
    data_1 = np.interp(energy, energy_1, signal_1)
    data_2 = np.interp(energy, energy_2, signal_2)
    return np.sqrt(np.sum(np.square(data_2 - data_1))) / len(energy)


def compare_spectra(energy_1, signal_1, energy_2, signal_2, initial_offset=0, initial_amplitude=1.0):
    """find best offset and amplitude"""

    def min_fun(vals):
        offset_val, amplitude_val = vals
        return spectra_difference(energy_1, signal_1, energy_2 + offset_val, signal_2 * amplitude_val)

    out = minimize(min_fun, np.array([initial_offset, initial_amplitude]), method='Powell')  # 'Nelder-Mead', 'Powell'

    offset, amplitude = out.x
    print('Minimise completed')
    print('Results:')
    print(out.message)
    print('Minimum = ', out.fun)
    print(f"  Offset = {offset:.2f} eV")
    print(f"  Amplitude = {amplitude:.3g}")
    return offset, amplitude
