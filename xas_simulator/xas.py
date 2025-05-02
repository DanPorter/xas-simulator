import os
import tempfile
import numpy as np
from scipy.integrate import romb
from nexpy.gui.datadialogs import GridParameters, NXDialog
from nexpy.gui.widgets import NXMessageBox
from nexpy.gui.utils import report_error, keep_data, display_message
from nexusformat.nexus import NeXusError, NXdata, NXfield, NXentry, NXroot

from xas_simulator.params_short import xray_data, parameters
from xas_simulator.simulation import DEFAULT_QUANTY_PATH, TMPDIR
from xas_simulator.XMCD_src2 import XAS_Lua
from .info_dialogue import INFODialog


def show_dialog():
    try:
        dialog = XASDialog()
        dialog.show()
    except NeXusError as error:
        report_error("XAS simulator setup", error)


class XASDialog(NXDialog):

    def __init__(self, parent=None):

        super().__init__(parent)

        # self.select_entry()
        self.parameters = GridParameters()
        self.parameters.add('Qtypath', DEFAULT_QUANTY_PATH, 'Quanty path')
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
        self.set_layout(self.parameters.grid(),
                        self.action_buttons(('Plot', self.plot_data)),
                        self.action_buttons(('Compute', self.compute)))
        self.set_title('XAS setup')

    def plot_data(self):
        # self.convert_QE().plot()
        pass

    def create_table(self, data):
        html_content = "<table> <tr>"
        nrow = len(data)
        ncol = len(data[0])
        # first rows are headers
        for icol in range(ncol):
            rtext = f'<th>{data[0][icol]}</th>'
            html_content += rtext
        html_content += '</tr>'
        # next rows are data:
        for irow in range(1, nrow):
            html_content += '<tr>'
            for icol in range(ncol):
                html_content += f'<td>{data[irow][icol]}</td>'
            html_content += '</tr>'
        html_content += '</table>'
        # print(html_content)
        return html_content

    def compute(self):

        qtypath = self.parameters['Qtypath'].value
        ion = self.parameters['ion'].value

        ch = int(self.parameters['charge'].value)
        sym = self.parameters['sym'].value
        label = ion + '_XAS'
        beta = self.parameters['beta'].value
        Dq = self.parameters['10Dq'].value
        mag_field = [self.parameters['Bx'].value,
                     self.parameters['By'].value,
                     self.parameters['Bz'].value]
        exchange_field = [self.parameters['Hx'].value,
                          self.parameters['Hy'].value,
                          self.parameters['Hz'].value]
        temperature = self.parameters['T'].value

        # Check ion
        if ion not in parameters or ion not in xray_data['elements']:
            message = f"Ion '{ion}' not available. Available ions are:\n"
            message += ', '.join(parameters)
            display_message('XAS Setup', message)
            return

        params2 = {
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
        # Check charge
        ch_str = f"{abs(ch)}+" if ch > 0 else f"{abs(ch)}-"
        if ch_str not in xray_data['elements'][ion]['charges']:
            message = f"Ionic charge: '{ion}{ch_str}' is not available.\nAvailable charges for {ion} are:\n"
            message += ','.join(xray_data['elements'][ion]['charges'].keys())
            display_message('XAS Setup', message)
            return

        # Check sym
        if sym not in xray_data['elements'][ion]['charges'][ch_str]['symmetries']:
            message = f"Symmetry: '{sym}' is not available for {ion}{ch_str}.\nAvailable symmetries for {ion}{ch_str} are:\n"
            message += ','.join(xray_data['elements'][ion]['charges'][ch_str]['symmetries'].keys())
            display_message('XAS Setup', message)
            return

        input = XAS_Lua(
            ion=ion,
            symm=sym,
            charge=ch_str,
            params=params2,
            params_json=xray_data,
            output_path=TMPDIR,
            quanty_path=qtypath,
        )
        input.write_header()
        input.H_init()
        input.setH_terms()
        input.set_electrons()
        input.define_atomic_term()
        input.define_crystal_field_term()
        input.define_external_field_term()
        input.setTemperature()
        input.setRestrictions()
        input.set_iterative_solver()
        input.set_spectra_functions()
        input.define_transitions([0, 0, 1], [0, 1, 0], [1, 0, 0])
        input.set_spectra_lists()
        # make and use a tmp directory
        input.calculate_and_save_spectra()
        input.run()
        input.treat_output()

        # Post-processing and display

        xt = np.loadtxt(os.path.join(input.path, label + '_iso.spec'), skiprows=5)
        mcd = np.loadtxt(os.path.join(input.path, label + '_cd.spec'), skiprows=5)
        xl = np.loadtxt(os.path.join(input.path, label + '_l.spec'), skiprows=5)
        xr = np.loadtxt(os.path.join(input.path, label + '_r.spec'), skiprows=5)

        xz = 3 * xt - (xl + xr)
        mcd2 = xr.copy()
        mcd2[:, 2] = xr[:, 2] - xl[:, 2]
        npts = np.shape(xz)[0]
        print(sum(abs(mcd2 - mcd)))
        # TOTAL spectra
        xas = xt.copy()
        xas[:, 2] = xz[:, 2] + xl[:, 2] + xr[:, 2]

        xas0 = xz.copy()
        xas0[:, 2] = (xl[:, 2] + xr[:, 2]) / 2 + xl[:, 2] + xr[:, 2]

        dx = xz.copy()
        dx[:, 2] = (xl[:, 2] + xr[:, 2]) - 2 * xz[:, 2]
        nh = 10 - (parameters[ion]['Nelec'] - (int(ch) - 2))
        # print(nh)
        trapz = True
        if trapz:
            tot = np.trapz(xas[:, 2], xas[:, 0])
            tot0 = np.trapz(xas0[:, 2], xas0[:, 0])
            dx0 = np.trapz(dx[:, 2], dx[:, 0])
        else:
            tot = romb(xas[:, 2], dx=float(xas[1, 0] - xas[0, 0]))
            tot0 = romb(xas0[:, 2], dx=float(xas0[1, 0] - xas0[0, 0]))
            dx0 = romb(dx[:, 2], dx=float(dx[1, 0] - dx[0, 0]))

        deltaXas = dx0 / tot

        if trapz:

            lz = 2 * nh * np.trapz(mcd2[:, 2], mcd2[:, 0]) / tot
            szef = 3 / 2 * nh * (
                        np.trapz(mcd2[0:npts // 2, 2], mcd2[0:npts // 2, 0]) - 2 * np.trapz(mcd2[npts // 2:, 2],
                                                                                            mcd2[npts // 2:, 0])) / tot
            lz0 = 2 * nh * np.trapz(mcd2[:, 2], mcd2[:, 0]) / tot0
            szef0 = 3 / 2 * nh * (
                        np.trapz(mcd2[0:npts // 2, 2], mcd2[0:npts // 2, 0]) - 2 * np.trapz(mcd2[npts // 2:, 2],
                                                                                            mcd2[npts // 2:, 0])) / tot0

        else:
            # print(len(mcd2[npts//2:,2]), len(mcd2[0:npts//2+1]))
            mydelta = mcd2[1, 0] - mcd2[0, 0]
            lz = 2 * nh * romb(mcd2[:, 2], mydelta) / tot
            L3 = romb(mcd2[0:npts // 2 + 1, 2], mydelta)
            L2 = romb(mcd2[npts // 2:, 2], mydelta)
            szef = 3 / 2 * nh * (L3 - 2 * L2) / tot
            lz0 = 2 * nh * romb(mcd2[:, 2], mydelta) / tot0
            szef0 = 3 / 2 * nh * (L3 - 2 * L2) / tot0

        Lz_t = float(input.outdic['L_k'])
        Sz_t = float(input.outdic['S_k'])
        Tz_t = float(input.outdic['T_k'])
        Seff_t = float(input.outdic['S_k']) + float(input.outdic['T_k'])
        #
        table1 = [[r'L$_z$', r'S$_{eff}$', r'S$_{z}$', r'T$_{z}$'],
                  [Lz_t, Seff_t, Sz_t, Tz_t]]

        # Create tables data with subscript HTML formatting
        table1_data = [['L<sub>z</sub>', 'S<sub>eff</sub>', 'S<sub>z</sub>', 'T<sub>z</sub>'],
                       [Lz_t, Seff_t, Sz_t, Tz_t]]
        table2_data = [['sL<sub>z</sub>', 'sS<sub>eff</sub>'],
                       [np.round(lz, decimals=4), np.round(szef, decimals=4)]]
        table3_data = [['s<sub>0</sub>L<sub>z</sub>', 's<sub>0</sub>S<sub>eff</sub>'],
                       [np.round(lz0, decimals=4), np.round(szef0, decimals=4)]]
        table4_data = [['ΔL<sub>z</sub> (%)', 'ΔS<sub>eff</sub> (%)', 'ΔXAS (%)', 'Δ<sub>0</sub>L<sub>z</sub> (%)',
                        'Δ<sub>0</sub>S<sub>eff</sub> (%)'],
                       [np.round(100 * (abs(Lz_t) - abs(lz)) / Lz_t, decimals=4),
                        np.round(100 * (abs(Seff_t) - abs(szef)) / Seff_t, decimals=4),
                        np.round(deltaXas * 100, decimals=4),
                        np.round(100 * (abs(Lz_t) - abs(lz0)) / Lz_t, decimals=4),
                        np.round(100 * (abs(Seff_t) - abs(szef0)) / Seff_t, decimals=4)]]

        html_content = self.create_table(table1_data)
        html_content += self.create_table(table2_data)
        html_content += self.create_table(table3_data)
        html_content += self.create_table(table4_data)
        resbox = NXMessageBox('XAS sum rules', html_content)
        resbox.setStyleSheet("QScrollArea{min-width:600 px; min-height: 400px}")
        resdiag = NXDialog()
        resdiag.set_layout(resbox)
        resdiag.set_title('Sum rules')
        resdiag.show()

        # for plotting
        Ene = NXfield(xt[:, 0], name='Energy')
        XAStot = NXfield(xt[:, 2], name='XAS total')
        XMCD = NXfield(mcd2[:, 2], name='XMCD')

        result1 = NXdata(XMCD, Ene)
        result2 = NXdata(XAStot, Ene)
        result1.plot('g-', label='XMCD')

        keep_data(result1)
        keep_data(result2)
