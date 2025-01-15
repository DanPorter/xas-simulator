import os
import tempfile

from xas_simulator.params_short import xray_data, parameters
from xas_simulator.XMCD_src2 import XAS_Lua

# DEFAULT_QUANTY_PATH = '/scratch/grp66007/software/xmcd_beamline_simulator/quanty_lin/Quanty'
# DEFAULT_QUANTY_PATH = r"C:\Users\grp66007\Documents\quanty\quanty_win\QuantyWin64.exe"
DEFAULT_QUANTY_PATH = 'Quanty'


# Find writable directory
TMPDIR = tempfile.gettempdir()
if not os.access(TMPDIR, os.W_OK):
    TMPDIR = os.path.abspath(os.path.join(os.path.dirname(__file__), '..'))
    if not os.access(TMPDIR, os.W_OK):
        TMPDIR = os.path.expanduser('~')
print('Writable TEMPDIR = %s' % TMPDIR)


def create_simulation(ion: str, ch: float, beta: float, dq: float, mag_field: tuple[float, float, float],
                      exchange_field: tuple[float, float, float], temperature: float,
                      quanty_path=DEFAULT_QUANTY_PATH) -> XAS_Lua:
    """
    Generate parameters for Quanty Simulation
    """

    # Check ion
    if ion not in parameters or ion not in xray_data['elements']:
        message = f"Ion '{ion}' not available. Available ions are:\n"
        message += ', '.join(parameters)
        raise Exception(message)

    # Check charge
    ch_str = f"{abs(ch)}+" if ch > 0 else f"{abs(ch)}-"
    if ch_str not in xray_data['elements'][ion]['charges']:
        message = f"Ionic charge: '{ion}{ch_str}' is not available.\nAvailable charges for {ion} are:\n"
        message += ','.join(xray_data['elements'][ion]['charges'].keys())
        raise Exception(message)

    # build parameters
    calculation_parameters = {
        'Nelec': parameters[ion]['Nelec'],
        'H_atomic': 1,
        'H_crystal_field': 1,
        'H_3d_ligands_hybridization_lmct': 0,
        'H_3d_ligands_hybridization_mlct': 0,
        'H_magnetic_field': 1,
        'H_exchange_field': 1,
        'Bx_i': mag_field[0],
        'By_i': mag_field[1],
        'Bz_i': mag_field[2],
        'Bx_f': mag_field[0],
        'By_f': mag_field[1],
        'Bz_f': mag_field[2],
        'Hx_i': exchange_field[0],
        'Hy_i': exchange_field[1],
        'Hz_i': exchange_field[2],
        'Hx_f': exchange_field[0],
        'Hy_f': exchange_field[1],
        'Hz_f': exchange_field[2],
        'T': temperature,
    }
    simulation = XAS_Lua(
        ion=ion,
        symm='Oh',
        charge=ch_str,
        params=calculation_parameters,
        params_json=xray_data,
        output_path=TMPDIR,
        quanty_path=quanty_path,
    )
    return simulation

