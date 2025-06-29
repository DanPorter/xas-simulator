import os
import subprocess


class XAS_Lua:
    """
    This class generates a lua file to be used as input by Quanty

    """

    def __init__(self, ion, symm, charge, params, params_json, output_path, quanty_path, beta=None):
        self.ion = ion
        fname = ion + '_XAS.lua'
        self.charge = charge
        self.path = output_path
        self.filename = os.path.join(output_path, fname)
        self.label = fname
        self.params = params
        self.xdat = params_json
        self.symm = symm
        self.Qty_path = quanty_path
        self.iondata = params_json['elements'][ion]['charges'][charge]['symmetries'][symm]['experiments']['XAS'][
            'edges']
        if beta is None:
            self.beta = {
                'F2dd_i': 0.8,
                'F2dd_f': 0.8,
                'F4dd_i': 0.8,
                'F4dd_f': 0.8,
                'zeta_3d': 1,
                'Xzeta_3d': 1,
                'zeta_2p': 1,
                'F2pd_f': 0.8,
                'G1pd_f': 0.8,
                'G3pd_f': 0.8
            }
        self.result = None
        self.outdic = None
        

    def write_header(self):
        """Write a header to the Lua file."""
        with open(self.filename, 'w') as f:
            f.write('-- This is an auto-generated Lua file.\n\n')

    def H_init(self):
        """Initialize the Hamiltonians in the Lua file."""
        with open(self.filename, 'a') as f:
            f.write(70 * '-' + '\n-- Initialize the Hamiltonians.\n' + 70 * '-' + '\n')
            f.write('H_i = 0\n')
            f.write('H_f = 0\n\n')

    def setH_terms(self):
        
        """Toggle the Hamiltonian terms based on provided parameters."""
        with open(self.filename, 'a') as f:
            f.write(70 * '-' + '\n-- Toggle the Hamiltonian terms.\n' + 70 * '-' + '\n')

            # Get values from the self.params dictionary; if not present, then, it's 0
            tog = ["H_atomic", "H_crystal_field", "H_3d_ligands_hybridization_lmct",
                   "H_3d_ligands_hybridization_mlct", "H_magnetic_field", "H_exchange_field"]
            for toggle in tog:
                dummy = self.params[toggle]
                if dummy is None:
                    self.params.update({toggle:'0'})
                    
            f.write(f'H_atomic = {self.params["H_atomic"]}\n')
            f.write(f'H_crystal_field = {self.params["H_crystal_field"]}\n')
            f.write(f'H_3d_ligands_hybridization_lmct = {self.params["H_3d_ligands_hybridization_lmct"]}\n')
            f.write(f'H_3d_ligands_hybridization_mlct = {self.params["H_3d_ligands_hybridization_mlct"]}\n')
            f.write(f'H_magnetic_field = {self.params["H_magnetic_field"]}\n')
            f.write(f'H_exchange_field = {self.params["H_exchange_field"]}\n\n')

    def set_electrons(self):
        """Set basic electronic information in the Lua file."""
        n3d = self.params['Nelec']
        with open(self.filename, 'a') as f:
            f.write(70 * '-' + '\n-- Basic information about electrons\n' + 70 * '-' + '\n')
            f.write('NBosons = 0\n')
            f.write('NFermions = 16\n\n')
            f.write('NElectrons_2p = 6\n')
            f.write('NElectrons_3d =%2d\n\n' % (n3d))
            f.write('IndexDn_2p = {0, 2, 4}\n')
            f.write('IndexUp_2p = {1, 3, 5}\n')
            f.write('IndexDn_3d = {6, 8, 10, 12, 14}\n')
            f.write('IndexUp_3d = {7, 9, 11, 13, 15}\n\n')

            if (self.params["H_3d_ligands_hybridization_lmct"]):
                f.write('NFermions = 26 \n')
                f.write('NElectrons_L1 = 10 \n')
                f.write('IndexDn_L1 = {16, 18, 20, 22, 24} \n')
                f.write('IndexUp_L1 = {17, 19, 21, 23, 25} \n')
            if(self.params["H_3d_ligands_hybridization_mlct"]):
                f.write('NFermions = 26 \n')
                f.write('NElectrons_L2 = 0 \n')
                f.write('IndexDn_L2 = {16, 18, 20, 22, 24} \n')
                f.write('IndexUp_L2 = {17, 19, 21, 23, 25} \n')
            

    def define_atomic_term(self):
        """
        Set the atomic part of the Hamiltonian
        """

        def get_atomic_params(ion, charge, xdat, Nelec):
            if str(int(charge[0]) - 2) != 0:
                conf = '3d' + str(Nelec - int(charge[0]) + 2)
                conf_xas = '2p5,3d' + str(Nelec - int(charge[0]) + 2 + 1)
            else:
                conf = '3d' + str(Nelec)
                conf_xas = '2p5,3d' + str(Nelec + 1)

            print(conf, conf_xas)
            indict = xdat['elements'][ion]['charges'][charge]['configurations'][conf]['terms']['Atomic']['parameters'][
                'variable']
            findict = \
                xdat['elements'][ion]['charges'][charge]['configurations'][conf_xas]['terms']['Atomic']['parameters'][
                    'variable']
            return [indict, findict]

        ini, fin = get_atomic_params(self.ion, self.charge, self.xdat, self.params['Nelec'])

        with open(self.filename, 'a') as f:
            f.write(70 * '-' + '\n-- Define the atomic term.\n' + 70 * '-' + '\n')

            # Define number operators for 2p and 3d orbitals
            f.write(
                "N_2p = NewOperator('Number', NFermions, IndexUp_2p, IndexUp_2p, {1, 1, 1})\n" +
                "+ NewOperator('Number', NFermions, IndexDn_2p, IndexDn_2p, {1, 1, 1})\n")
            f.write(
                "N_3d = NewOperator('Number', NFermions, IndexUp_3d, IndexUp_3d, {1, 1, 1, 1, 1})\n"+
                "+  NewOperator('Number', NFermions, IndexDn_3d, IndexDn_3d, {1, 1, 1, 1, 1})\n")

            # Check condition for H_atomic
            if self.params.get('H_atomic', 1) == 1:
                f.write("    F0_3d_3d = NewOperator('U', NFermions, IndexUp_3d, IndexDn_3d, {1, 0, 0})\n")
                f.write("    F2_3d_3d = NewOperator('U', NFermions, IndexUp_3d, IndexDn_3d, {0, 1, 0})\n")
                f.write("    F4_3d_3d = NewOperator('U', NFermions, IndexUp_3d, IndexDn_3d, {0, 0, 1})\n")

                f.write(
                    "    F0_2p_3d = NewOperator('U', NFermions, IndexUp_2p, IndexDn_2p, IndexUp_3d, IndexDn_3d, {1, 0}, {0, 0})\n")
                f.write(
                    "    F2_2p_3d = NewOperator('U', NFermions, IndexUp_2p, IndexDn_2p, IndexUp_3d, IndexDn_3d, {0, 1}, {0, 0})\n")
                f.write(
                    "    G1_2p_3d = NewOperator('U', NFermions, IndexUp_2p, IndexDn_2p, IndexUp_3d, IndexDn_3d, {0, 0}, {1, 0})\n")
                f.write(
                    "    G3_2p_3d = NewOperator('U', NFermions, IndexUp_2p, IndexDn_2p, IndexUp_3d, IndexDn_3d, {0, 0}, {0, 1})\n")
                # initial
                f.write(f"    U_3d_3d_i = {ini['U(3d,3d)']}\n")
                f.write(f"    F2_3d_3d_i = {ini['F2(3d,3d)']} * {self.beta['F2dd_i']}\n")
                f.write(f"    F4_3d_3d_i = {ini['F4(3d,3d)']} * {self.beta['F4dd_i']}\n")
                f.write(f"    F0_3d_3d_i = U_3d_3d_i + 2 / 63 * F2_3d_3d_i + 2 / 63 * F4_3d_3d_i\n")
                # final
                f.write(f"    U_3d_3d_f = {fin['U(3d,3d)']}\n")
                f.write(f"    F2_3d_3d_f = {fin['F2(3d,3d)']} * {self.beta['F2dd_f']}\n")
                f.write(f"    F4_3d_3d_f = {fin['F4(3d,3d)']} * {self.beta['F4dd_f']}\n")
                f.write(f"    F0_3d_3d_f = U_3d_3d_f + 2 / 63 * F2_3d_3d_f + 2 / 63 * F4_3d_3d_f\n")
                f.write(f"    U_2p_3d_f =  {fin['U(2p,3d)']}\n")
                f.write(f"    F2_2p_3d_f = {fin['F2(2p,3d)']} * {self.beta['F2pd_f']}\n")
                f.write(f"    G1_2p_3d_f = {fin['G1(2p,3d)']} * {self.beta['G1pd_f']}\n")
                f.write(f"    G3_2p_3d_f = {fin['G3(2p,3d)']} * {self.beta['G3pd_f']}\n")
                f.write(f"    F0_2p_3d_f = U_2p_3d_f + 1 / 15 * G1_2p_3d_f + 3 / 70 * G3_2p_3d_f\n")

                f.write(
                    "    H_i = H_i + Chop( F0_3d_3d_i * F0_3d_3d + F2_3d_3d_i * F2_3d_3d + F4_3d_3d_i * F4_3d_3d)\n")
                f.write(
                    "    H_f = H_f + Chop( F0_3d_3d_f * F0_3d_3d + F2_3d_3d_f * F2_3d_3d + F4_3d_3d_f * F4_3d_3d "+
                    " + F0_2p_3d_f * F0_2p_3d + F2_2p_3d_f * F2_2p_3d + G1_2p_3d_f * G1_2p_3d + G3_2p_3d_f * G3_2p_3d)\n")

                f.write("    ldots_3d = NewOperator('ldots', NFermions, IndexUp_3d, IndexDn_3d)\n")
                f.write("    ldots_2p = NewOperator('ldots', NFermions, IndexUp_2p, IndexDn_2p)\n")

                f.write(f"    zeta_3d_i = {ini['ζ(3d)'] * self.beta['zeta_3d']}\n")
                f.write(f"    zeta_3d_f = {fin['ζ(3d)'] * self.beta['Xzeta_3d']}\n")
                f.write(f"    zeta_2p_f = {fin['ζ(2p)'] * self.beta['zeta_2p']}\n")

                f.write("    H_i = H_i + Chop( zeta_3d_i * ldots_3d)\n")
                f.write("    H_f = H_f + Chop( zeta_3d_f * ldots_3d + zeta_2p_f * ldots_2p)\n\n")



        
    def define_crystal_field_term(self):
        """
        Set the crystal field part of the Hamiltonian
        """
        Nelec = self.params['Nelec']
        if str(int(self.charge[0]) - 2) != 0:
            conf = '3d' + str(Nelec - int(self.charge[0]) + 2)
            conf_xas = '2p5,3d' + str(Nelec - int(self.charge[0]) + 2 + 1)
        else:
            conf = '3d' + str(Nelec)
            conf_xas = '2p5,3d' + str(Nelec + 1)

        match self.symm:
            case 'Oh':
                self.define_Oh_crystal_field_term(conf, conf_xas)
            case 'D3h':
                self.define_D3h_crystal_field_term(conf,conf_xas)
            case 'D4h':
                self.define_D4h_crystal_field_term(conf, conf_xas)
            case 'Td':
                self.define_Td_crystal_field_term(conf,conf_xas)
            case 'C3v':
                self.define_C3v_crystal_field_term(conf,conf_xas)
            case _:
                self.define_Oh_crystal_field_term(conf,conf_xas)

    def define_Oh_crystal_field_term(self, conf, conf_xas):
        """
        Crystal field with Oh symmetry
        """
        # Attempt 1: Check for separate initial and final values
        tendq_i = self.params.get('10Dq_i')
        tendq_f = self.params.get('10Dq_f')

        if tendq_i is None or tendq_f is None:
            # Attempt 2: Check for a single 10Dq value
            single_tendq = self.params.get('10Dq')
            if single_tendq is not None:
                tendq_i = single_tendq
                tendq_f = single_tendq
            else:
                # Attempt 3: Fallback to the deeply nested structure
                try:
                    tendq_i = self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf]['terms']['Crystal Field']['symmetries']['Oh']['parameters']['variable']['10Dq(3d)']
                    tendq_f = self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf_xas]['terms']['Crystal Field']['symmetries']['Oh']['parameters']['variable']['10Dq(3d)']
                except KeyError:
                    print("Error: 10Dq parameter not found in any of the expected locations.")

                    
        with open(self.filename, 'a') as f:
            f.write(70 * '-' + '\n-- Define the crystal field term.\n' + 70 * '-' + '\n')
            if self.params.get('H_crystal_field', 1) == 1:
                f.write('Akm = {{4, 0, 2.1}, {4, -4, 1.5 * sqrt(0.7)}, {4, 4, 1.5 * sqrt(0.7)}}\n')
                f.write("tenDq_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, Akm)\n")
                f.write(f"tenDq_3d_i = {tendq_i} \n")
                f.write(f"tenDq_3d_f = {tendq_f} \n")

                f.write("H_i = H_i + Chop(tenDq_3d_i * tenDq_3d)\n")
                f.write("H_f = H_f + Chop(tenDq_3d_f * tenDq_3d)\n\n")
        if (self.params["H_3d_ligands_hybridization_lmct"]):
            self.define_Oh_crystal_field_lmct(conf, conf_xas)
        if (self.params["H_3d_ligands_hybridization_mlct"]):
            self.define_Oh_crystal_field_mlct(conf, conf_xas)


    def define_Oh_crystal_field_lmct(self, conf, conf_xas):
        """
        Ligand field for Oh symmetry
        Delta_3d_L1
        U_3d_3d
        E_L1
        tenDq_L1
        Veg_3d_L1
        Vt2g_3d_L1
        """
        terms = ['Delta_L1',  'tenDq_L1', 'Veg_L1', 'Vt2g_L1']

        for it in terms:
            dummy_i = self.params.get(it+'_i')
            dummy_f = self.params.get(it+'_f')
            
            if dummy_i is None or dummy_f is None:
                dummy = self.params.get(it)
                if dummy is None:
                    self.params.update({it:'0'})
                else:
                    self.params.update({it+'_i':dummy})
                    self.params.update({it+'_f':dummy})


        with open(self.filename, 'a') as f:
            f.write('    N_L1 = NewOperator("Number", NFermions, IndexUp_L1, IndexUp_L1, {1, 1, 1, 1, 1})\n')
            f.write('         + NewOperator("Number", NFermions, IndexDn_L1, IndexDn_L1, {1, 1, 1, 1, 1})\n')
            f.write('\n')
            f.write(f'    Delta_3d_L1_i = {self.params["Delta_L1_i"]}\n')
            f.write('    E_3d_i = (10 * Delta_3d_L1_i - NElectrons_3d * (19 + NElectrons_3d) * U_3d_3d_i / 2) / (10 + NElectrons_3d)\n')
            f.write('    E_L1_i = NElectrons_3d * ((1 + NElectrons_3d) * U_3d_3d_i / 2 - Delta_3d_L1_i) / (10 + NElectrons_3d)\n')
            f.write('\n')
            f.write(f'    Delta_3d_L1_f = {self.params["Delta_L1_f"]} \n')
            f.write('    E_3d_f = (10 * Delta_3d_L1_f - NElectrons_3d * (31 + NElectrons_3d) * U_3d_3d_f / 2 - 90 * U_2p_3d_f) / (16 + NElectrons_3d)\n')
            f.write('    E_2p_f = (10 * Delta_3d_L1_f + (1 + NElectrons_3d) * (NElectrons_3d * U_3d_3d_f / 2 - (10 + NElectrons_3d) * U_2p_3d_f)) / (16 + NElectrons_3d)\n')
            f.write('    E_L1_f = ((1 + NElectrons_3d) * (NElectrons_3d * U_3d_3d_f / 2 + 6 * U_2p_3d_f) - (6 + NElectrons_3d) * Delta_3d_L1_f) / (16 + NElectrons_3d)\n')
            f.write('\n')
            f.write('    H_i = H_i + Chop(\n')
            f.write('          E_3d_i * N_3d\n')
            f.write('        + E_L1_i * N_L1)\n')
            f.write('\n')
            f.write('    H_f = H_f + Chop(\n')
            f.write('          E_3d_f * N_3d\n')
            f.write('        + E_2p_f * N_2p\n')
            f.write('        + E_L1_f * N_L1)\n')
            f.write('\n')
            f.write('    tenDq_L1 = NewOperator("CF", NFermions, IndexUp_L1, IndexDn_L1, PotentialExpandedOnClm("Oh", 2, {0.6, -0.4}))\n')
            f.write('\n')
            f.write('    Veg_3d_L1 = NewOperator("CF", NFermions, IndexUp_L1, IndexDn_L1, IndexUp_3d, IndexDn_3d, PotentialExpandedOnClm("Oh", 2, {1, 0}))\n')
            f.write('              + NewOperator("CF", NFermions, IndexUp_3d, IndexDn_3d, IndexUp_L1, IndexDn_L1, PotentialExpandedOnClm("Oh", 2, {1, 0}))\n')
            f.write('\n')
            f.write('    Vt2g_3d_L1 = NewOperator("CF", NFermions, IndexUp_L1, IndexDn_L1, IndexUp_3d, IndexDn_3d, PotentialExpandedOnClm("Oh", 2, {0, 1}))\n')
            f.write('               + NewOperator("CF", NFermions, IndexUp_3d, IndexDn_3d, IndexUp_L1, IndexDn_L1, PotentialExpandedOnClm("Oh", 2, {0, 1}))\n')
            f.write('\n')
            f.write(f'    tenDq_L1_i = {self.params["tenDq_L1_i"]} \n')
            f.write(f'    Veg_3d_L1_i = {self.params["Veg_L1_i"]} \n')
            f.write(f'    Vt2g_3d_L1_i = {self.params["Vt2g_L1_i"]} \n')
            f.write('\n')
            f.write(f'    tenDq_L1_f = {self.params["tenDq_L1_f"]} \n')
            f.write(f'    Veg_3d_L1_f = {self.params["Veg_L1_f"]} \n')
            f.write(f'    Vt2g_3d_L1_f ={self.params["Vt2g_L1_f"]} \n')
            f.write('\n')
            f.write('    H_i = H_i + Chop(\n')
            f.write('          tenDq_L1_i * tenDq_L1\n')
            f.write('        + Veg_3d_L1_i * Veg_3d_L1\n')
            f.write('        + Vt2g_3d_L1_i * Vt2g_3d_L1)\n')
            f.write('\n')
            f.write('    H_f = H_f + Chop(\n')
            f.write('          tenDq_L1_f * tenDq_L1\n')
            f.write('        + Veg_3d_L1_f * Veg_3d_L1\n')
            f.write('        + Vt2g_3d_L1_f * Vt2g_3d_L1)')

    def define_Oh_crystal_field_mlct(self, conf, conf_xas):
        """
        Ligand field for Oh symmetry
        Delta_3d_L2
        U_3d_3d
        E_L2
        tenDq_L2
        Veg_3d_L2
        Vt2g_3d_L2
        """
        terms = ['Delta_L2',  'tenDq_L2', 'Veg_L2', 'Vt2g_L2']

        for it in terms:
            dummy_i = self.params.get(it+'_i')
            dummy_f = self.params.get(it+'_f')
            
            if dummy_i is None or dummy_f is None:
                dummy = self.params.get(it)
                if dummy is None:
                    self.params.update({it:'0'})
                else:
                    self.params.update({it+'_i':dummy})
                    self.params.update({it+'_f':dummy})


        with open(self.filename, 'a') as f:
            f.write('    N_L2 = NewOperator("Number", NFermions, IndexUp_L2, IndexUp_L2, {1, 1, 1, 1, 1})\n')
            f.write('         + NewOperator("Number", NFermions, IndexDn_L2, IndexDn_L2, {1, 1, 1, 1, 1})\n')
            f.write('\n')
            f.write(f'    Delta_3d_L2_i = {self.params["Delta_L2_i"]}\n')
            f.write('    E_3d_i = U_3d_3d_i * (-NElectrons_3d + 1) / 2\n')
            f.write('    E_L2_i = Delta_3d_L2_i + U_3d_3d_i * NElectrons_3d / 2 - U_3d_3d_i / 2 \n')
            f.write('\n')
            f.write(f'    Delta_3d_L2_f = {self.params["Delta_L2_f"]} \n')
            f.write('    E_3d_f = -(U_3d_3d_f * NElectrons_3d^2 + 11 * U_3d_3d_f * NElectrons_3d + 60 * U_2p_3d_f) / (2 * NElectrons_3d + 12) \n')
            f.write('    E_2p_f = NElectrons_3d * (U_3d_3d_f * NElectrons_3d + U_3d_3d_f - 2 * U_2p_3d_f * NElectrons_3d - 2 * U_2p_3d_f) / (2 * (NElectrons_3d + 6)) \n')
            f.write('    E_L2_f = (2 * Delta_3d_L2_f * NElectrons_3d + 12 * Delta_3d_L2_f + U_3d_3d_f * NElectrons_3d^2 - U_3d_3d_f * NElectrons_3d - 12 * U_3d_3d_f + 12 * U_2p_3d_f * NElectrons_3d + 12 * U_2p_3d_f) / (2 * (NElectrons_3d + 6)) \n')
            f.write('\n')
            f.write('    H_i = H_i + Chop(\n')
            f.write('          E_3d_i * N_3d\n')
            f.write('        + E_L2_i * N_L2)\n')
            f.write('\n')
            f.write('    H_f = H_f + Chop(\n')
            f.write('          E_3d_f * N_3d\n')
            f.write('        + E_2p_f * N_2p\n')
            f.write('        + E_L2_f * N_L2)\n')
            f.write('\n')
            f.write('    tenDq_L2 = NewOperator("CF", NFermions, IndexUp_L2, IndexDn_L2, PotentialExpandedOnClm("Oh", 2, {0.6, -0.4}))\n')
            f.write('\n')
            f.write('    Veg_3d_L2 = NewOperator("CF", NFermions, IndexUp_L2, IndexDn_L2, IndexUp_3d, IndexDn_3d, PotentialExpandedOnClm("Oh", 2, {1, 0}))\n')
            f.write('              + NewOperator("CF", NFermions, IndexUp_3d, IndexDn_3d, IndexUp_L2, IndexDn_L2, PotentialExpandedOnClm("Oh", 2, {1, 0}))\n')
            f.write('\n')
            f.write('    Vt2g_3d_L2 = NewOperator("CF", NFermions, IndexUp_L2, IndexDn_L2, IndexUp_3d, IndexDn_3d, PotentialExpandedOnClm("Oh", 2, {0, 1}))\n')
            f.write('               + NewOperator("CF", NFermions, IndexUp_3d, IndexDn_3d, IndexUp_L2, IndexDn_L2, PotentialExpandedOnClm("Oh", 2, {0, 1}))\n')
            f.write('\n')
            f.write(f'    tenDq_L2_i = {self.params["tenDq_L2_i"]} \n')
            f.write(f'    Veg_3d_L2_i = {self.params["Veg_L2_i"]} \n')
            f.write(f'    Vt2g_3d_L2_i = {self.params["Vt2g_L2_i"]} \n')
            f.write('\n')
            f.write(f'    tenDq_L2_f = {self.params["tenDq_L2_f"]} \n')
            f.write(f'    Veg_3d_L2_f = {self.params["Veg_L2_f"]} \n')
            f.write(f'    Vt2g_3d_L2_f ={self.params["Vt2g_L2_f"]} \n')
            f.write('\n')
            f.write('    H_i = H_i + Chop(\n')
            f.write('          tenDq_L2_i * tenDq_L2\n')
            f.write('        + Veg_3d_L2_i * Veg_3d_L2\n')
            f.write('        + Vt2g_3d_L2_i * Vt2g_3d_L2)\n')
            f.write('\n')
            f.write('    H_f = H_f + Chop(\n')
            f.write('          tenDq_L2_f * tenDq_L2\n')
            f.write('        + Veg_3d_L2_f * Veg_3d_L2\n')
            f.write('        + Vt2g_3d_L2_f * Vt2g_3d_L2)')

    def define_D3h_crystal_field_term(self, conf, conf_xas):
        """
        Crystal field with D3h symmetry
        """
        # Attempt 1: Check for separate initial and final values
        Dmu_i = self.params.get('Dmu_i')
        Dmu_f = self.params.get('Dmu_f')

        if Dmu_i is None or Dmu_f is None:
            # Attempt 2: Check for a single Dmu value
            single_Dmu = self.params.get('Dmu')
            if single_Dmu is not None:
                Dmu_i = single_Dmu
                Dmu_f = single_Dmu
            else:
                # Attempt 3: Fallback to the json data
                try:
                    Dmu_i =  self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf]['terms']['Crystal Field']['symmetries']['D3h']['parameters']['variable']['Dμ(3d)']
                    Dmu_f =  self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf_xas]['terms']['Crystal Field']['symmetries']['D3h']['parameters']['variable']['Dμ(3d)']
                except KeyError:
                    print("Error: Dmu parameter not found in any of the expected locations.")

        Dnu_i = self.params.get('Dnu_i')
        Dnu_f = self.params.get('Dnu_f')

        if Dnu_i is None or Dnu_f is None:
            # Attempt 2: Check for a single Dnu value
            single_Dnu = self.params.get('Dnu')
            if single_Dnu is not None:
                Dnu_i = single_Dnu
                Dnu_f = single_Dnu
            else:
                # Attempt 3: Fallback to the json data
                try:
                    Dnu_i =  self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf]['terms']['Crystal Field']['symmetries']['D3h']['parameters']['variable']['Dν(3d)']
                    Dnu_f =  self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf_xas]['terms']['Crystal Field']['symmetries']['D3h']['parameters']['variable']['Dν(3d)']
                except KeyError:
                    print("Error: Dnu parameter not found in any of the expected locations.")

        with open(self.filename, 'a') as f:
            f.write(70 * '-' + '\n-- Define the crystal field term.\n' + 70 * '-' + '\n')
            if self.params.get('H_crystal_field', 1) == 1:
                f.write('Akm = {{2, 0, -7}}\n')
                f.write("Dmu_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, Akm)\n")
                f.write('Akm = {{4, 0, -21}}\n')
                f.write("Dnu_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, Akm)\n")
                f.write(f"Dmu_3d_i = {Dmu_i} \n")
                f.write(f"Dmu_3d_f = {Dmu_f} \n")
                f.write(f"Dnu_3d_i = {Dnu_i} \n")
                f.write(f"Dnu_3d_f = {Dnu_f} \n")

                f.write("H_i = H_i + Chop(Dmu_3d_i * Dmu_3d + Dnu_3d_i * Dnu_3d)\n")
                f.write("H_f = H_f + Chop(Dmu_3d_f * Dmu_3d + Dnu_3d_f * Dnu_3d)\n\n")
        

        
    def define_D4h_crystal_field_term(self, conf, conf_xas):
        """
        Crystal field with D4h symmetry
        'Dq(3d)', 'Ds(3d)', 'Dt(3d)'
        """
       # Attempt 1: Check for separate initial and final values
        Dq_i = self.params.get('Dq_i')
        Dq_f = self.params.get('Dq_f')

        if Dq_i is None or Dq_f is None:
            # Attempt 2: Check for a single Dq value
            single_Dq = self.params.get('Dq')
            if single_Dq is not None:
                Dq_i = single_Dq
                Dq_f = single_Dq
            else:
                # Attempt 3: Fallback to the json data
                try:
                    Dq_i =  self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf]['terms']['Crystal Field']['symmetries']['D4h']['parameters']['variable']['Dq(3d)']
                    Dq_f =  self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf_xas]['terms']['Crystal Field']['symmetries']['D4h']['parameters']['variable']['Dq(3d)']

                except KeyError:
                    print("Error: Dq parameter not found in any of the expected locations.")

       # Attempt 1: Check for separate initial and final values
        Ds_i = self.params.get('Ds_i')
        Ds_f = self.params.get('Ds_f')

        if Ds_i is None or Ds_f is None:
            # Attempt 2: Check for a single Ds value
            single_Ds = self.params.get('Ds')
            if single_Ds is not None:
                Ds_i = single_Ds
                Ds_f = single_Ds
            else:
                # Attempt 3: Fallback to the json data
                try:
                    Ds_i =  self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf]['terms']['Crystal Field']['symmetries']['D4h']['parameters']['variable']['Ds(3d)']
                    Ds_f =  self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf_xas]['terms']['Crystal Field']['symmetries']['D4h']['parameters']['variable']['Ds(3d)']
                except KeyError:
                    print("Error: Ds parameter not found in any of the expected locations.")

       # Attempt 1: Check for separate initial and final values
        Dt_i = self.params.get('Dt_i')
        Dt_f = self.params.get('Dt_f')

        if Dt_i is None or Dt_f is None:
            # Attempt 2: Check for a single Dt value
            single_Dt = self.params.get('Dt')
            if single_Dt is not None:
                Dt_i = single_Dt
                Dt_f = single_Dt
            else:
                # Attempt 3: Fallback to the json data
                try:
                    Dt_i =  self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf]['terms']['Crystal Field']['symmetries']['D4h']['parameters']['variable']['Dt(3d)']
                    Dt_f =  self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf_xas]['terms']['Crystal Field']['symmetries']['D4h']['parameters']['variable']['Dt(3d)']
                except KeyError:
                    print("Error: Dt parameter not found in any of the expected locations.")                    
                    


        with open(self.filename, 'a') as f:
            f.write(70 * '-' + '\n-- Define the crystal field term.\n' + 70 * '-' + '\n')
            if self.params.get('H_crystal_field', 1) == 1:
                f.write('Akm = {{4, 0, 21}, {4, -4, 1.5 * sqrt(70)}, {4, 4, 1.5 * sqrt(70)}}\n')
                f.write("Dq_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, Akm)\n")
                f.write('Akm = {{2, 0, -7}}\n')
                f.write("Ds_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, Akm)\n")
                f.write('Akm = {{4, 0, -21}}\n')
                f.write("Dt_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, Akm)\n")
                
                f.write(f"Dq_3d_i = {Dq_i} \n")
                f.write(f"Dq_3d_f = {Dq_f} \n")
                f.write(f"Ds_3d_i = {Ds_i} \n")
                f.write(f"Ds_3d_f = {Ds_f} \n")
                f.write(f"Dt_3d_i = {Dt_i} \n")
                f.write(f"Dt_3d_f = {Dt_f} \n")

                f.write("H_i = H_i + Chop(Dq_3d_i * Dq_3d + Ds_3d_i * Ds_3d + Dt_3d_i * Dt_3d)\n")
                f.write("H_f = H_f + Chop(Dq_3d_f * Dq_3d + Ds_3d_f * Ds_3d + Dt_3d_f * Dt_3d)\n")

        if (self.params["H_3d_ligands_hybridization_lmct"]):
            self.define_D4h_crystal_field_lmct(conf, conf_xas)
        if (self.params["H_3d_ligands_hybridization_mlct"]):
            self.define_D4h_crystal_field_mlct(conf, conf_xas)


        
    def define_D4h_crystal_field_lmct(self, conf, conf_xas):
        """
        Delta_L1
        Va1g_L1
        Vb1g_L1
        Vb2g_L1
        Veg_L1
        Dq_L1
        Ds_L1
        Dt_L1
        """
        terms = ['Delta_L1','Va1g_L1','Vb1g_L1','Vb2g_L1','Veg_L1','Dq_L1','Ds_L1','Dt_L1']

        for it in terms:
            dummy_i = self.params.get(it+'_i')
            dummy_f = self.params.get(it+'_f')
            
            if dummy_i is None or dummy_f is None:
                dummy = self.params.get(it)
                if dummy is None:
                    self.params.update({it:'0'})
                else:
                    self.params.update({it+'_i':dummy})
                    self.params.update({it+'_f':dummy})
        with open(self.filename, 'a') as f:
            f.write('N_L1 = NewOperator("Number", NFermions, IndexUp_L1, IndexUp_L1, {1, 1, 1, 1, 1}) \n')
            f.write('     + NewOperator("Number", NFermions, IndexDn_L1, IndexDn_L1, {1, 1, 1, 1, 1}) \n')
            f.write(' \n')
            f.write(f'Delta_3d_L1_i = {self.params["Delta_L1_i"]} \n')
            f.write('E_3d_i = (10 * Delta_3d_L1_i - NElectrons_3d * (19 + NElectrons_3d) * U_3d_3d_i / 2) / (10 + NElectrons_3d) \n')
            f.write('E_L1_i = NElectrons_3d * ((1 + NElectrons_3d) * U_3d_3d_i / 2 - Delta_3d_L1_i) / (10 + NElectrons_3d) \n')
            f.write(' \n')
            f.write(f'Delta_3d_L1_f = {self.params["Delta_L1_f"]} \n')
            f.write('E_3d_f = (10 * Delta_3d_L1_f - NElectrons_3d * (31 + NElectrons_3d) * U_3d_3d_f / 2 - 90 * U_2p_3d_f) / (16 + NElectrons_3d) \n')
            f.write('E_2p_f = (10 * Delta_3d_L1_f + (1 + NElectrons_3d) * (NElectrons_3d * U_3d_3d_f / 2 - (10 + NElectrons_3d) * U_2p_3d_f)) / (16 + NElectrons_3d) \n')
            f.write('E_L1_f = ((1 + NElectrons_3d) * (NElectrons_3d * U_3d_3d_f / 2 + 6 * U_2p_3d_f) - (6 + NElectrons_3d) * Delta_3d_L1_f) / (16 + NElectrons_3d) \n')
            f.write(' \n')
            f.write('H_i = H_i + Chop( \n')
            f.write('      E_3d_i * N_3d \n')
            f.write('    + E_L1_i * N_L1) \n')
            f.write(' \n')
            f.write('H_f = H_f + Chop( \n')
            f.write('      E_3d_f * N_3d \n')
            f.write('    + E_2p_f * N_2p \n')
            f.write('    + E_L1_f * N_L1) \n')
            f.write(' \n')
            f.write('Dq_L1 = NewOperator("CF", NFermions, IndexUp_L1, IndexDn_L1, PotentialExpandedOnClm("D4h", 2, { 6,  6, -4, -4})) \n')
            f.write('Ds_L1 = NewOperator("CF", NFermions, IndexUp_L1, IndexDn_L1, PotentialExpandedOnClm("D4h", 2, {-2,  2,  2, -1})) \n')
            f.write('Dt_L1 = NewOperator("CF", NFermions, IndexUp_L1, IndexDn_L1, PotentialExpandedOnClm("D4h", 2, {-6, -1, -1,  4})) \n')
            f.write(' \n')
            f.write('Va1g_3d_L1 = NewOperator("CF", NFermions, IndexUp_L1, IndexDn_L1, IndexUp_3d, IndexDn_3d, PotentialExpandedOnClm("D4h", 2, {1, 0, 0, 0})) \n')
            f.write('           + NewOperator("CF", NFermions, IndexUp_3d, IndexDn_3d, IndexUp_L1, IndexDn_L1, PotentialExpandedOnClm("D4h", 2, {1, 0, 0, 0})) \n')
            f.write(' \n')
            f.write('Vb1g_3d_L1 = NewOperator("CF", NFermions, IndexUp_L1, IndexDn_L1, IndexUp_3d, IndexDn_3d, PotentialExpandedOnClm("D4h", 2, {0, 1, 0, 0})) \n')
            f.write('           + NewOperator("CF", NFermions, IndexUp_3d, IndexDn_3d, IndexUp_L1, IndexDn_L1, PotentialExpandedOnClm("D4h", 2, {0, 1, 0, 0})) \n')
            f.write(' \n')
            f.write('Vb2g_3d_L1 = NewOperator("CF", NFermions, IndexUp_L1, IndexDn_L1, IndexUp_3d, IndexDn_3d, PotentialExpandedOnClm("D4h", 2, {0, 0, 1, 0})) \n')
            f.write('           + NewOperator("CF", NFermions, IndexUp_3d, IndexDn_3d, IndexUp_L1, IndexDn_L1, PotentialExpandedOnClm("D4h", 2, {0, 0, 1, 0})) \n')
            f.write(' \n')
            f.write('Veg_3d_L1 = NewOperator("CF", NFermions, IndexUp_L1, IndexDn_L1, IndexUp_3d, IndexDn_3d, PotentialExpandedOnClm("D4h", 2, {0, 0, 0, 1})) \n')
            f.write('          + NewOperator("CF", NFermions, IndexUp_3d, IndexDn_3d, IndexUp_L1, IndexDn_L1, PotentialExpandedOnClm("D4h", 2, {0, 0, 0, 1})) \n')
            f.write(' \n')
            f.write(f'Dq_L1_i = {self.params["Dq_L1_i"]} \n')
            f.write(f'Ds_L1_i = {self.params["Ds_L1_i"]} \n')
            f.write(f'Dt_L1_i = {self.params["Dt_L1_i"]} \n')
            f.write(f'Va1g_3d_L1_i = {self.params["Va1g_L1_i"]} \n')
            f.write(f'Vb1g_3d_L1_i = {self.params["Vb1g_L1_i"]} \n')
            f.write(f'Vb2g_3d_L1_i = {self.params["Vb2g_L1_i"]} \n')
            f.write(f'Veg_3d_L1_i =  {self.params["Veg_L1_i"]}\n')
            f.write(' \n')
            f.write(f'Dq_L1_f = {self.params["Dq_L1_f"]} \n')
            f.write(f'Ds_L1_f = {self.params["Ds_L1_f"]} \n')
            f.write(f'Dt_L1_f = {self.params["Dt_L1_f"]} \n')
            f.write(f'Va1g_3d_L1_f = {self.params["Va1g_L1_f"]} \n')
            f.write(f'Vb1g_3d_L1_f = {self.params["Vb1g_L1_f"]} \n')
            f.write(f'Vb2g_3d_L1_f = {self.params["Vb2g_L1_f"]} \n')
            f.write(f'Veg_3d_L1_f = {self.params["Veg_L1_f"]} \n')
            f.write(' \n')
            f.write('H_i = H_i + Chop( \n')
            f.write('      Dq_L1_i * Dq_L1 \n')
            f.write('    + Ds_L1_i * Ds_L1 \n')
            f.write('    + Dt_L1_i * Dt_L1 \n')
            f.write('    + Va1g_3d_L1_i * Va1g_3d_L1 \n')
            f.write('    + Vb1g_3d_L1_i * Vb1g_3d_L1 \n')
            f.write('    + Vb2g_3d_L1_i * Vb2g_3d_L1 \n')
            f.write('    + Veg_3d_L1_i  * Veg_3d_L1) \n')
            f.write(' \n')
            f.write('H_f = H_f + Chop( \n')
            f.write('      Dq_L1_f * Dq_L1 \n')
            f.write('    + Ds_L1_f * Ds_L1 \n')
            f.write('    + Dt_L1_f * Dt_L1 \n')
            f.write('    + Va1g_3d_L1_f * Va1g_3d_L1 \n')
            f.write('    + Vb1g_3d_L1_f * Vb1g_3d_L1 \n')
            f.write('    + Vb2g_3d_L1_f * Vb2g_3d_L1 \n')
            f.write('    + Veg_3d_L1_f  * Veg_3d_L1) \n')

    def define_D4h_crystal_field_mlct(self, conf, conf_xas):
        """
        Delta_L2
        Va1g_L2
        Vb1g_L2
        Vb2g_L2
        Veg_L2
        Dq_L2
        Ds_L2
        Dt_L2
        """
        terms = ['Delta_L2','Va1g_L2','Vb1g_L2','Vb2g_L2','Veg_L2','Dq_L2','Ds_L2','Dt_L2']

        for it in terms:
            dummy_i = self.params.get(it+'_i')
            dummy_f = self.params.get(it+'_f')
            
            if dummy_i is None or dummy_f is None:
                dummy = self.params.get(it)
                if dummy is None:
                    self.params.update({it:'0'})
                else:
                    self.params.update({it+'_i':dummy})
                    self.params.update({it+'_f':dummy})
        with open(self.filename, 'a') as f:
            f.write('N_L2 = NewOperator("Number", NFermions, IndexUp_L2, IndexUp_L2, {1, 1, 1, 1, 1}) \n')
            f.write('     + NewOperator("Number", NFermions, IndexDn_L2, IndexDn_L2, {1, 1, 1, 1, 1}) \n')
            f.write(' \n')
            f.write(f'Delta_3d_L2_i = {self.params["Delta_L2_i"]} \n')
            f.write(f'Delta_3d_L2_f = {self.params["Delta_L2_f"]} \n')
            f.write('E_3d_i = U_3d_3d_i * (-NElectrons_3d + 1) / 2 \n')
            f.write('E_L2_i = Delta_3d_L2_i + U_3d_3d_i * NElectrons_3d / 2 - U_3d_3d_i / 2 \n')
            f.write(' \n')
            f.write('E_3d_f = -(U_3d_3d_f * NElectrons_3d^2 + 11 * U_3d_3d_f * NElectrons_3d + 60 * U_2p_3d_f) / (2 * NElectrons_3d + 12) \n')
            f.write('E_2p_f = NElectrons_3d * (U_3d_3d_f * NElectrons_3d + U_3d_3d_f - 2 * U_2p_3d_f * NElectrons_3d - 2 * U_2p_3d_f) / (2 * (NElectrons_3d + 6)) \n')
            f.write('E_L2_f = (2 * Delta_3d_L2_f * NElectrons_3d + 12 * Delta_3d_L2_f + U_3d_3d_f * NElectrons_3d^2 - U_3d_3d_f * NElectrons_3d - 12 * U_3d_3d_f + 12 * U_2p_3d_f * NElectrons_3d + 12 * U_2p_3d_f) / (2 * (NElectrons_3d + 6)) \n')
            f.write(' \n')

            f.write('H_i = H_i + Chop( \n')
            f.write('      E_3d_i * N_3d \n')
            f.write('    + E_L2_i * N_L2) \n')
            f.write(' \n')
            f.write('H_f = H_f + Chop( \n')
            f.write('      E_3d_f * N_3d \n')
            f.write('    + E_2p_f * N_2p \n')
            f.write('    + E_L2_f * N_L2) \n')
            f.write(' \n')
            f.write('Dq_L2 = NewOperator("CF", NFermions, IndexUp_L2, IndexDn_L2, PotentialExpandedOnClm("D4h", 2, { 6,  6, -4, -4})) \n')
            f.write('Ds_L2 = NewOperator("CF", NFermions, IndexUp_L2, IndexDn_L2, PotentialExpandedOnClm("D4h", 2, {-2,  2,  2, -1})) \n')
            f.write('Dt_L2 = NewOperator("CF", NFermions, IndexUp_L2, IndexDn_L2, PotentialExpandedOnClm("D4h", 2, {-6, -1, -1,  4})) \n')
            f.write(' \n')
            f.write('Va1g_3d_L2 = NewOperator("CF", NFermions, IndexUp_L2, IndexDn_L2, IndexUp_3d, IndexDn_3d, PotentialExpandedOnClm("D4h", 2, {1, 0, 0, 0})) \n')
            f.write('           + NewOperator("CF", NFermions, IndexUp_3d, IndexDn_3d, IndexUp_L2, IndexDn_L2, PotentialExpandedOnClm("D4h", 2, {1, 0, 0, 0})) \n')
            f.write(' \n')
            f.write('Vb1g_3d_L2 = NewOperator("CF", NFermions, IndexUp_L2, IndexDn_L2, IndexUp_3d, IndexDn_3d, PotentialExpandedOnClm("D4h", 2, {0, 1, 0, 0})) \n')
            f.write('           + NewOperator("CF", NFermions, IndexUp_3d, IndexDn_3d, IndexUp_L2, IndexDn_L2, PotentialExpandedOnClm("D4h", 2, {0, 1, 0, 0})) \n')
            f.write(' \n')
            f.write('Vb2g_3d_L2 = NewOperator("CF", NFermions, IndexUp_L2, IndexDn_L2, IndexUp_3d, IndexDn_3d, PotentialExpandedOnClm("D4h", 2, {0, 0, 1, 0})) \n')
            f.write('           + NewOperator("CF", NFermions, IndexUp_3d, IndexDn_3d, IndexUp_L2, IndexDn_L2, PotentialExpandedOnClm("D4h", 2, {0, 0, 1, 0})) \n')
            f.write(' \n')
            f.write('Veg_3d_L2 = NewOperator("CF", NFermions, IndexUp_L2, IndexDn_L2, IndexUp_3d, IndexDn_3d, PotentialExpandedOnClm("D4h", 2, {0, 0, 0, 1})) \n')
            f.write('          + NewOperator("CF", NFermions, IndexUp_3d, IndexDn_3d, IndexUp_L2, IndexDn_L2, PotentialExpandedOnClm("D4h", 2, {0, 0, 0, 1})) \n')
            f.write(' \n')
            f.write(f'Dq_L2_i = {self.params["Dq_L2_i"]} \n')
            f.write(f'Ds_L2_i = {self.params["Ds_L2_i"]} \n')
            f.write(f'Dt_L2_i = {self.params["Dt_L2_i"]} \n')
            f.write(f'Va1g_3d_L2_i = {self.params["Va1g_L2_i"]} \n')
            f.write(f'Vb1g_3d_L2_i = {self.params["Vb1g_L2_i"]} \n')
            f.write(f'Vb2g_3d_L2_i = {self.params["Vb2g_L2_i"]} \n')
            f.write(f'Veg_3d_L2_i =  {self.params["Veg_L2_i"]}\n')
            f.write(' \n')
            f.write(f'Dq_L2_f = {self.params["Dq_L2_f"]} \n')
            f.write(f'Ds_L2_f = {self.params["Ds_L2_f"]} \n')
            f.write(f'Dt_L2_f = {self.params["Dt_L2_f"]} \n')
            f.write(f'Va1g_3d_L2_f = {self.params["Va1g_L2_f"]} \n')
            f.write(f'Vb1g_3d_L2_f = {self.params["Vb1g_L2_f"]} \n')
            f.write(f'Vb2g_3d_L2_f = {self.params["Vb2g_L2_f"]} \n')
            f.write(f'Veg_3d_L2_f = {self.params["Veg_L2_f"]} \n')
            f.write(' \n')
            f.write('H_i = H_i + Chop( \n')
            f.write('      Dq_L2_i * Dq_L2 \n')
            f.write('    + Ds_L2_i * Ds_L2 \n')
            f.write('    + Dt_L2_i * Dt_L2 \n')
            f.write('    + Va1g_3d_L2_i * Va1g_3d_L2 \n')
            f.write('    + Vb1g_3d_L2_i * Vb1g_3d_L2 \n')
            f.write('    + Vb2g_3d_L2_i * Vb2g_3d_L2 \n')
            f.write('    + Veg_3d_L2_i  * Veg_3d_L2) \n')
            f.write(' \n')
            f.write('H_f = H_f + Chop( \n')
            f.write('      Dq_L2_f * Dq_L2 \n')
            f.write('    + Ds_L2_f * Ds_L2 \n')
            f.write('    + Dt_L2_f * Dt_L2 \n')
            f.write('    + Va1g_3d_L2_f * Va1g_3d_L2 \n')
            f.write('    + Vb1g_3d_L2_f * Vb1g_3d_L2 \n')
            f.write('    + Vb2g_3d_L2_f * Vb2g_3d_L2 \n')
            f.write('    + Veg_3d_L2_f  * Veg_3d_L2) \n')

        
    def define_Td_crystal_field_term(self, conf, conf_xas):
        """
        Crystal field with Td symmetry
        """

        # Attempt 1: Check for separate initial and final values
        tendq_i = self.params.get('Dq_i')
        tendq_f = self.params.get('Dq_f')

        if tendq_i is None or tendq_f is None:
            # Attempt 2: Check for a single 10Dq value
            single_tendq = self.params.get('Dq')
            if single_tendq is not None:
                tendq_i = single_tendq
                tendq_f = single_tendq
            else:
                # Attempt 3: Fallback to the deeply nested structure
                try:
                    tendq_i = self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf]['terms']['Crystal Field']['symmetries']['Td']['parameters']['variable']['10Dq(3d)']
                    tendq_f = self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf_xas]['terms']['Crystal Field']['symmetries']['Td']['parameters']['variable']['10Dq(3d)']
                except KeyError:
                    print("Error: 10Dq parameter not found in any of the expected locations.")

        with open(self.filename, 'a') as f:
            f.write(70 * '-' + '\n-- Define the crystal field term.\n' + 70 * '-' + '\n')
            if self.params.get('H_crystal_field', 1) == 1:
                f.write('Akm = {{4, 0, -2.1}, {4, -4, -1.5 * sqrt(0.7)}, {4, 4, -1.5 * sqrt(0.7)}}\n')
                f.write("tenDq_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, Akm)\n")
                f.write(f"tenDq_3d_i = {tendq_i} \n")
                f.write(f"tenDq_3d_f = {tendq_f} \n")

                f.write("H_i = H_i + Chop(tenDq_3d_i * tenDq_3d)\n")
                f.write("H_f = H_f + Chop(tenDq_3d_f * tenDq_3d)\n\n")

    def define_C3v_crystal_field_term(self, conf, conf_xas):
        """
        Crystal field with C3v symmetry
        'Dq(3d)', 'Dσ(3d)', 'Dτ(3d)'
        """
       # Attempt 1: Check for separate initial and final values
        Dq_i = self.params.get('Dq_i')
        Dq_f = self.params.get('Dq_f')

        if Dq_i is None or Dq_f is None:
            # Attempt 2: Check for a single Dq value
            single_Dq = self.params.get('Dq')
            if single_Dq is not None:
                Dq_i = single_Dq
                Dq_f = single_Dq
            else:
                # Attempt 3: Fallback to the json data
                try:
                    Dq_i =  self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf]['terms']['Crystal Field']['symmetries']['C3v']['parameters']['variable']['Dq(3d)']
                    Dq_f =  self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf_xas]['terms']['Crystal Field']['symmetries']['C3v']['parameters']['variable']['Dq(3d)']
                except KeyError:
                    print("Error: Dq parameter not found in any of the expected locations.")


       # Attempt 1: Check for separate initial and final values
        Dsigma_i = self.params.get('Dsigma_i')
        Dsigma_f = self.params.get('Dsigma_f')

        if Dsigma_i is None or Dsigma_f is None:
            # Attempt 2: Check for a single Dsigma value
            single_Dsigma = self.params.get('Dsigma')
            if single_Dsigma is not None:
                Dsigma_i = single_Dsigma
                Dsigma_f = single_Dsigma
            else:
                # Attempt 3: Fallback to the json data
                try:
                    Dsigma_i =  self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf]['terms']['Crystal Field']['symmetries']['C3v']['parameters']['variable']['Dσ(3d)']
                    Dsigma_f =  self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf_xas]['terms']['Crystal Field']['symmetries']['C3v']['parameters']['variable']['Dσ(3d)']
                except KeyError:
                    print("Error: Dsigma parameter not found in any of the expected locations.")
                    
       # Attempt 1: Check for separate initial and final values
        Dtau_i = self.params.get('Dtau_i')
        Dtau_f = self.params.get('Dtau_f')

        if Dtau_i is None or Dtau_f is None:
            # Attempt 2: Check for a single Dtau value
            single_Dtau = self.params.get('Dtau')
            if single_Dtau is not None:
                Dtau_i = single_Dtau
                Dtau_f = single_Dtau
            else:
                # Attempt 3: Fallback to the json data
                try:
                    Dtau_i =  self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf]['terms']['Crystal Field']['symmetries']['C3v']['parameters']['variable']['Dτ(3d)']
                    Dtau_f =  self.xdat['elements'][self.ion]['charges'][self.charge]\
                        ['configurations'][conf_xas]['terms']['Crystal Field']['symmetries']['C3v']['parameters']['variable']['Dτ(3d)']
                except KeyError:
                    print("Error: Dtau parameter not found in any of the expected locations.")

        with open(self.filename, 'a') as f:
            f.write(70 * '-' + '\n-- Define the crystal field term.\n' + 70 * '-' + '\n')
            if self.params.get('H_crystal_field', 1) == 1:
                f.write('Akm = {{4, 0, -14}, {4, 3, -2 * math.sqrt(70)}, {4, -3, 2 * math.sqrt(70)}}\n')
                f.write("Dq_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, Akm)\n")
                f.write('Akm = {{2, 0, -7}}\n')
                f.write("Dsigma_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, Akm)\n")
                f.write('Akm = {{4, 0, -21}}\n')
                f.write("Dtau_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, Akm)\n")
                
                f.write(f"Dq_3d_i = {Dq_i} \n")
                f.write(f"Dq_3d_f = {Dq_f} \n")
                f.write(f"Dsigma_3d_i = {Dsigma_i} \n")
                f.write(f"Dsigma_3d_f = {Dsigma_f} \n")
                f.write(f"Dtau_3d_i = {Dtau_i} \n")
                f.write(f"Dtau_3d_f = {Dtau_f} \n")

                f.write("H_i = H_i + Chop(Dq_3d_i * Dq_3d + Dsigma_3d_i * Dsigma_3d + Dtau_3d_i * Dtau_3d)\n")
                f.write("H_f = H_f + Chop(Dq_3d_f * Dq_3d + Dsigma_3d_f * Dsigma_3d + Dtau_3d_f * Dtau_3d)\n")



    def define_external_field_term(self):
        """
        Set the external (magnetic, exchange) field part of the Hamiltonian
        """
        with open(self.filename, 'a') as f:
            f.write(70 * '-' + '\n-- Define the magnetic and exchange field terms.\n' + 70 * '-' + '\n')

            f.write("Sx_3d = NewOperator('Sx', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Sy_3d = NewOperator('Sy', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Sz_3d = NewOperator('Sz', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Ssqr_3d = NewOperator('Ssqr', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Splus_3d = NewOperator('Splus', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Smin_3d = NewOperator('Smin', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Lx_3d = NewOperator('Lx', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Ly_3d = NewOperator('Ly', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Lz_3d = NewOperator('Lz', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Lsqr_3d = NewOperator('Lsqr', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Lplus_3d = NewOperator('Lplus', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Lmin_3d = NewOperator('Lmin', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Jx_3d = NewOperator('Jx', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Jy_3d = NewOperator('Jy', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Jz_3d = NewOperator('Jz', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Jsqr_3d = NewOperator('Jsqr', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Jplus_3d = NewOperator('Jplus', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Jmin_3d = NewOperator('Jmin', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Tx_3d = NewOperator('Tx', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Ty_3d = NewOperator('Ty', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Tz_3d = NewOperator('Tz', NFermions, IndexUp_3d, IndexDn_3d)\n")
            f.write("Sx = Sx_3d\n")
            f.write("Sy = Sy_3d\n")
            f.write("Sz = Sz_3d\n")
            f.write("Lx = Lx_3d\n")
            f.write("Ly = Ly_3d\n")
            f.write("Lz = Lz_3d\n")
            f.write("Jx = Jx_3d\n")
            f.write("Jy = Jy_3d\n")
            f.write("Jz = Jz_3d\n")
            f.write("Tx = Tx_3d\n")
            f.write("Ty = Ty_3d\n")
            f.write("Tz = Tz_3d\n")
            f.write("Ssqr = Sx * Sx + Sy * Sy + Sz * Sz\n")
            f.write("Lsqr = Lx * Lx + Ly * Ly + Lz * Lz\n")
            f.write("Jsqr = Jx * Jx + Jy * Jy + Jz * Jz\n")
            f.write(f"Bx_i = {self.params['Bx_i']} * EnergyUnits.Tesla.value\n")
            f.write(f"By_i = {self.params['By_i']} * EnergyUnits.Tesla.value\n")
            f.write(f"Bz_i = {self.params['Bz_i']} * EnergyUnits.Tesla.value\n")
            f.write(f"Bx_f = {self.params['Bx_f']} * EnergyUnits.Tesla.value\n")
            f.write(f"By_f = {self.params['By_f']} * EnergyUnits.Tesla.value\n")
            f.write(f"Bz_f = {self.params['Bz_f']} * EnergyUnits.Tesla.value\n")
            f.write(f"H_i = H_i + Chop( Bx_i * (2 * Sx + Lx) + By_i * (2 * Sy + Ly) + Bz_i * (2 * Sz + Lz))\n")
            f.write(f"H_f = H_f + Chop( Bx_f * (2 * Sx + Lx) + By_f * (2 * Sy + Ly) + Bz_f * (2 * Sz + Lz))\n")
            f.write(f"Hx_i = {self.params['Hx_i']}\n")
            f.write(f"Hy_i = {self.params['Hy_i']}\n")
            f.write(f"Hz_i = {self.params['Hz_i']}\n")
            f.write(f"Hx_f = {self.params['Hx_f']}\n")
            f.write(f"Hy_f = {self.params['Hy_f']}\n")
            f.write(f"Hz_f = {self.params['Hz_f']}\n")
            f.write(f"H_i = H_i + Chop( Hx_i * Sx + Hy_i * Sy + Hz_i * Sz)\n")
            f.write(f"H_f = H_f + Chop( Hx_f * Sx + Hy_f * Sy + Hz_f * Sz)\n\n")

    def setTemperature(self):
        """
        set temperature
        """
        with open(self.filename, 'a') as f:
            f.write(70 * '-' + '\n-- Define the temperature.\n' + 70 * '-' + '\n')
            f.write(f"T={self.params['T']} * EnergyUnits.Kelvin.value\n\n")

    def setRestrictions(self):
        """
        """
        nconfs = self.params.get('NConfigurations')
        if nconfs is None:
            nconfs = 2
        with open(self.filename, 'a') as f:
            f.write(70 * '-' + '\n-- Define the restrictions.\n' + 70 * '-' + '\n')
            f.write("InitialRestrictions = {NFermions, NBosons, {'111111 0000000000', NElectrons_2p, NElectrons_2p},\n")
            f.write("                               {'000000 1111111111', NElectrons_3d, NElectrons_3d}}\n")

            f.write(
                "FinalRestrictions = {NFermions, NBosons, {'111111 0000000000', NElectrons_2p - 1, NElectrons_2p - 1},\n")
            f.write("                             {'000000 1111111111', NElectrons_3d + 1, NElectrons_3d + 1}}\n\n")
            if self.params['H_3d_ligands_hybridization_lmct']:
                
                f.write(f"NConfigurations = {nconfs} \n")
                f.write(
                    "    InitialRestrictions = {NFermions, NBosons, {'111111 0000000000 0000000000', NElectrons_2p, NElectrons_2p},\n")
                f.write(
                    "                                               {'000000 1111111111 0000000000', NElectrons_3d, NElectrons_3d},\n")
                f.write(
                    "                                               {'000000 0000000000 1111111111', NElectrons_L1, NElectrons_L1}}\n\n")
                f.write(
                    "    FinalRestrictions = {NFermions, NBosons, {'111111 0000000000 0000000000', NElectrons_2p - 1, NElectrons_2p - 1},\n")
                f.write(
                    "                                             {'000000 1111111111 0000000000', NElectrons_3d + 1, NElectrons_3d + 1},\n")
                f.write(
                    "                                             {'000000 0000000000 1111111111', NElectrons_L1, NElectrons_L1}}\n\n")
                f.write(
                    "    CalculationRestrictions = {NFermions, NBosons, {'000000 0000000000 1111111111', NElectrons_L1 - (NConfigurations - 1), NElectrons_L1}}\n\n")
            if self.params['H_3d_ligands_hybridization_mlct']:
                f.write(
                    "    InitialRestrictions = {NFermions, NBosons, {'111111 0000000000 0000000000', NElectrons_2p, NElectrons_2p},\n")
                f.write(
                    "                                               {'000000 1111111111 0000000000', NElectrons_3d, NElectrons_3d},\n")
                f.write(
                    "                                               {'000000 0000000000 1111111111', NElectrons_L2, NElectrons_L2}}\n\n")
                f.write(
                    "    FinalRestrictions = {NFermions, NBosons, {'111111 0000000000 0000000000', NElectrons_2p - 1, NElectrons_2p - 1},\n")
                f.write(
                    "                                             {'000000 1111111111 0000000000', NElectrons_3d + 1, NElectrons_3d + 1},\n")
                f.write(
                    "                                             {'000000 0000000000 1111111111', NElectrons_L2, NElectrons_L2}}\n\n")
                f.write(
                    "    CalculationRestrictions = {NFermions, NBosons, {'000000 0000000000 1111111111', NElectrons_L2, NElectrons_L2 + (NConfigurations - 1)}}\n\n")

    def set_iterative_solver(self):
        """
        solve
        """

        with open(self.filename, 'a') as f:
            f.write(70 * '-' + '\n-- Converge the number of states\n' + 70 * '-' + '\n')
            f.write("epsilon = 1.19e-07\n")
            f.write("\n")
            f.write("NPsis = 45.0\n")
            f.write("NPsisAuto = 1\n")
            f.write("\n")
            f.write("dZ = {}\n")
            f.write("\n")
            f.write("if NPsisAuto == 1 and NPsis ~= 1 then\n")
            f.write("    NPsis = 4\n")
            f.write("    NPsisIncrement = 8\n")
            f.write("    NPsisIsConverged = false\n")
            f.write("\n")
            f.write("    while not NPsisIsConverged do\n")
            f.write("        if CalculationRestrictions == nil then\n")
            f.write("            Psis_i = Eigensystem(H_i, InitialRestrictions, NPsis)\n")
            f.write("        else\n")
            f.write(
                "            Psis_i = Eigensystem(H_i, InitialRestrictions, NPsis, {{'restrictions', CalculationRestrictions}})\n")
            f.write("        end\n")
            f.write("\n")
            f.write("        if not (type(Psis_i) == 'table') then\n")
            f.write("            Psis_i = {Psis_i}\n")
            f.write("        end\n")
            f.write("\n")
            f.write("        E_gs_i = Psis_i[1] * H_i * Psis_i[1]\n")
            f.write("\n")
            f.write("        Z = 0\n")
            f.write("\n")
            f.write("        for i, Psi in ipairs(Psis_i) do\n")
            f.write("            E = Psi * H_i * Psi\n")
            f.write("\n")
            f.write("            if math.abs(E - E_gs_i) < epsilon then\n")
            f.write("                dZ[i] = 1\n")
            f.write("            else\n")
            f.write("                dZ[i] = math.exp(-(E - E_gs_i) / T)\n")
            f.write("            end\n")
            f.write("\n")
            f.write("            Z = Z + dZ[i]\n")
            f.write("\n")
            f.write("            if (dZ[i] / Z) < math.sqrt(epsilon) then\n")
            f.write("                i = i - 1\n")
            f.write("                NPsisIsConverged = true\n")
            f.write("                NPsis = i\n")
            f.write("                Psis_i = {unpack(Psis_i, 1, i)}\n")
            f.write("                dZ = {unpack(dZ, 1, i)}\n")
            f.write("                break\n")
            f.write("            end\n")
            f.write("        end\n")
            f.write("\n")
            f.write("        if NPsisIsConverged then\n")
            f.write("            break\n")
            f.write("        else\n")
            f.write("            NPsis = NPsis + NPsisIncrement\n")
            f.write("        end\n")
            f.write("    end\n")
            f.write("else\n")
            f.write("    if CalculationRestrictions == nil then\n")
            f.write("        Psis_i = Eigensystem(H_i, InitialRestrictions, NPsis)\n")
            f.write("    else\n")
            f.write(
                "        Psis_i = Eigensystem(H_i, InitialRestrictions, NPsis, {{'restrictions', CalculationRestrictions}})\n")
            f.write("    end\n")
            f.write("\n")
            f.write("    if not (type(Psis_i) == 'table') then\n")
            f.write("        Psis_i = {Psis_i}\n")
            f.write("    end\n")
            f.write("        E_gs_i = Psis_i[1] * H_i * Psis_i[1]\n")
            f.write("\n")
            f.write("    Z = 0\n")
            f.write("\n")
            f.write("    for i, Psi in ipairs(Psis_i) do\n")
            f.write("        E = Psi * H_i * Psi\n")
            f.write("\n")
            f.write("        if math.abs(E - E_gs_i) < epsilon then\n")
            f.write("            dZ[i] = 1\n")
            f.write("        else\n")
            f.write("            dZ[i] = math.exp(-(E - E_gs_i) / T)\n")
            f.write("        end\n")
            f.write("\n")
            f.write("        Z = Z + dZ[i]\n")
            f.write("    end\n")
            f.write("end\n")
            f.write("\n")
            f.write("-- Normalize dZ to unity.\n")
            f.write("for i in ipairs(dZ) do\n")
            f.write("    dZ[i] = dZ[i] / Z\n")
            f.write("end\n")

    def set_spectra_functions(self):
        """
        """
        iondata = \
            self.xdat['elements'][self.ion]['charges'][self.charge]['symmetries'][self.symm]['experiments']['XAS'][
                'edges']
        Edge = iondata['L2,3 (2p)']['axes'][0][4]
        Gmin = iondata['L2,3 (2p)']['axes'][0][5][0]
        Gmax = iondata['L2,3 (2p)']['axes'][0][5][1]
        Gamma = iondata['L2,3 (2p)']['axes'][0][6]
        Egamma1 = Edge + 10
        BaseName = self.path.replace('\\', '/') + '/' + self.ion + '_XAS'

        with open(self.filename, 'a') as f:
            f.write(70 * '-' + '\n-- Helper functions for spectra calculations\n' + 70 * '-' + '\n')
            f.write("function ValueInTable(value, table)\n")
            f.write("    -- Check if a value is in a table.\n")
            f.write("    for k, v in ipairs(table) do\n")
            f.write("        if value == v then\n")
            f.write("            return true\n")
            f.write("        end\n")
            f.write("    end\n")
            f.write("    return false\n")
            f.write("end\n")
            f.write("function GetSpectrum(G, T, Psis, indices, dZSpectra)\n")
            f.write("    -- Extract the spectra corresponding to the operators identified\n")
            f.write("    -- using the indices argument. The returned spectrum is a weighted\n")
            f.write("    -- sum, where the weights are the Boltzmann probabilities.\n")
            f.write("    if not (type(indices) == 'table') then\n")
            f.write("        indices = {" + " indices}\n")
            f.write("    end\n")
            f.write("    c = 1\n")
            f.write("    dZSpectrum = {}\n")
            f.write("    for i in ipairs(T) do\n")
            f.write("        for k in ipairs(Psis) do\n")
            f.write("            if ValueInTable(i, indices) then\n")
            f.write("                table.insert(dZSpectrum, dZSpectra[c])\n")
            f.write("            else\n")
            f.write("                table.insert(dZSpectrum, 0)\n")
            f.write("            end\n")
            f.write("            c = c + 1\n")
            f.write("        end\n")
            f.write("    end\n")
            f.write("    return Spectra.Sum(G, dZSpectrum)\n")
            f.write("end\n")
            f.write("function SaveSpectrum(G, suffix)\n")
            f.write("    -- Scale, broaden, and save the spectrum to disk.\n")
            f.write("    G = -1 / math.pi * G\n")
            f.write(f"    Gmin1 = {Gmin} - {Gamma}\n")
            f.write(f"    Gmax1 = {Gmax} - {Gamma}\n")
            f.write(f"    Egamma1 = ({Egamma1} - {Edge}) + DeltaE\n")
            f.write("    G.Broaden(0, {{Emin, Gmin1}, {Egamma1, Gmin1}, {Egamma1, Gmax1}, {Emax, Gmax1}})\n")
            f.write("    G.Print({{'file', '" + f"{BaseName}" + "_' .. suffix .. '.spec'}})\n")
            f.write("end\n\n")

    def define_transitions(self, kvec, eps11, eps12):
        """
        """
        # polarization information:

        with open(self.filename, 'a') as f:
            f.write(70 * '-' + '\n-- Define transitions \n' + 70 * '-' + '\n')
            f.write("t = math.sqrt(1/2)\n")
            f.write(
                "Tx_2p_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, IndexUp_2p, IndexDn_2p, {{1, -1, t    }, {1, 1, -t    }})\n")
            f.write(
                "Ty_2p_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, IndexUp_2p, IndexDn_2p, {{1, -1, t * I}, {1, 1,  t * I}})\n")
            f.write(
                "Tz_2p_3d = NewOperator('CF', NFermions, IndexUp_3d, IndexDn_3d, IndexUp_2p, IndexDn_2p, {{1,  0, 1    }                })\n")
            f.write("k = { %2d, %2d, %2d }\n" % (kvec[0], kvec[1], kvec[2]))
            f.write("ev = { %2d, %2d, %2d }\n" % (eps11[0], eps11[1], eps11[2]))
            f.write("eh = { %2d, %2d, %2d }\n" % (eps12[0], eps12[1], eps12[2]))
            # f.write(f"ev = {eps11}\n")
            # f.write(f"eh = {eps12}\n")
            f.write("-- Calculate the right and left polarization vectors.\n")
            f.write("er = {t * (eh[1] - I * ev[1]),\n")
            f.write("      t * (eh[2] - I * ev[2]),\n")
            f.write("      t * (eh[3] - I * ev[3])}\n")
            f.write("el = {-t * (eh[1] + I * ev[1]),\n")
            f.write("      -t * (eh[2] + I * ev[2]),\n")
            f.write("      -t * (eh[3] + I * ev[3])}\n")
            f.write("function CalculateT(e)\n")
            f.write("    -- Calculate the transition operator for arbitrary polarization.\n")
            f.write("    T = e[1] * Tx_2p_3d + e[2] * Ty_2p_3d + e[3] * Tz_2p_3d\n")
            f.write("    return Chop(T)\n")
            f.write("end\n")
            f.write("Tv_2p_3d = CalculateT(ev)\n")
            f.write("Th_2p_3d = CalculateT(eh)\n")
            f.write("Tr_2p_3d = CalculateT(er)\n")
            f.write("Tl_2p_3d = CalculateT(el)\n")
            f.write("Tk_2p_3d = CalculateT(k)\n\n")

    def set_spectra_lists(self):
        """
        """
        with open(self.filename, 'a') as f:
            f.write(70 * '-' + '\n-- List of spectra \n' + 70 * '-' + '\n')
            f.write("spectra = {'Isotropic','Circular Dichroism', 'Linear Dichroism'}\n")
            f.write("-- Create two lists, one with the operators and the second with\n")
            f.write("-- the indices of the operators required to calculate a given\n")
            f.write("-- spectrum.\n")
            f.write("T_2p_3d = {}\n")
            f.write("indices_2p_3d = {}\n")
            f.write("c = 1\n")
            f.write("spectrum = 'Isotropic'\n")
            f.write("if ValueInTable(spectrum, spectra) then\n")
            f.write("    indices_2p_3d[spectrum] = {}\n")
            f.write("    for j, operator in ipairs({Tr_2p_3d, Tl_2p_3d, Tk_2p_3d}) do\n")
            f.write("        table.insert(T_2p_3d, operator)\n")
            f.write("        table.insert(indices_2p_3d[spectrum], c)\n")
            f.write("        c = c + 1\n")
            f.write("    end\n")
            f.write("end\n")
            f.write("spectrum = 'Circular Dichroism'\n")
            f.write("if ValueInTable(spectrum, spectra) then\n")
            f.write("    indices_2p_3d[spectrum] = {}\n")
            f.write("    if ValueInTable('Isotropic', spectra) then\n")
            f.write("        table.insert(indices_2p_3d[spectrum], 1)\n")
            f.write("        table.insert(indices_2p_3d[spectrum], 2)\n")
            f.write("    else\n")
            f.write("        for j, operator in ipairs({Tr_2p_3d, Tl_2p_3d}) do\n")
            f.write("            table.insert(T_2p_3d, operator)\n")
            f.write("            table.insert(indices_2p_3d[spectrum], c)\n")
            f.write("            c = c + 1\n")
            f.write("        end\n")
            f.write("    end\n")
            f.write("end\n")
            f.write("spectrum = 'Linear Dichroism'\n")
            f.write("if ValueInTable(spectrum, spectra) then\n")
            f.write("    indices_2p_3d[spectrum] = {}\n")
            f.write("    for j, operator in ipairs({Tv_2p_3d, Th_2p_3d}) do\n")
            f.write("        table.insert(T_2p_3d, operator)\n")
            f.write("        table.insert(indices_2p_3d[spectrum], c)\n")
            f.write("        c = c + 1\n")
            f.write("    end\n")
            f.write("end\n\n")

    def calculate_and_save_spectra(self):
        """

        """
        element_data = self.xdat['elements'][self.ion]['charges'][self.charge]
        iondata = element_data['symmetries'][self.symm]['experiments']['XAS']['edges']
        Edge = iondata['L2,3 (2p)']['axes'][0][4]
        Gmin = iondata['L2,3 (2p)']['axes'][0][5][0]
        Gmax = iondata['L2,3 (2p)']['axes'][0][5][1]
        Gamma = iondata['L2,3 (2p)']['axes'][0][6]
        Egamma1 = Edge + 10
        Emin1 = Edge - 10
        Emax1 = Edge + 30
        NE1 = 2048
        DenseBorder = 2048

        with open(self.filename, 'a') as f:
            f.write(70 * '-' + '\n-- Calculate and save spectra \n' + 70 * '-' + '\n')
            f.write("Sk = Chop(k[1] * Sx + k[2] * Sy + k[3] * Sz)\n")
            f.write("Lk = Chop(k[1] * Lx + k[2] * Ly + k[3] * Lz)\n")
            f.write("Jk = Chop(k[1] * Jx + k[2] * Jy + k[3] * Jz)\n")
            f.write("Tk = Chop(k[1] * Tx + k[2] * Ty + k[3] * Tz)\n")
            f.write("Operators = {H_i, Ssqr, Lsqr, Jsqr, Sk, Lk, Jk, Tk, ldots_3d, N_2p, N_3d, 'dZ'}\n")
            f.write(r"header = 'Analysis of the initial Hamiltonian:\n'")
            f.write("\n")
            f.write(
                r"header = header .. '=================================================================================================================================\n'")
            f.write("\n")
            f.write(
                r"header = header .. 'State         <E>     <S^2>     <L^2>     <J^2>      <Sk>      <Lk>      <Jk>      <Tk>     <l.s>    <N_2p>    <N_3d>          dZ\n'")
            f.write("\n")
            f.write(
                r"header = header .. '=================================================================================================================================\n'")
            f.write("\n")
            f.write(
                r"footer = '=================================================================================================================================\n'")
            f.write("\n")
            f.write("if H_3d_ligands_hybridization_lmct == 1 then\n")
            f.write("    Operators = {H_i, Ssqr, Lsqr, Jsqr, Sk, Lk, Jk, Tk, ldots_3d, N_2p, N_3d, N_L1, 'dZ'}\n")
            f.write(r"    header = 'Analysis of the initial Hamiltonian:\n'")
            f.write("\n")
            f.write(
                r"    header = header .. '===========================================================================================================================================\n'")
            f.write("\n")
            f.write(
                r"    header = header .. 'State         <E>     <S^2>     <L^2>     <J^2>      <Sk>      <Lk>      <Jk>      <Tk>     <l.s>    <N_2p>    <N_3d>    <N_L1>          dZ\n'")
            f.write("\n")
            f.write(
                r"    header = header .. '===========================================================================================================================================\n'")
            f.write("\n")
            f.write(
                r"    footer = '=========================================================================================================================================== \n \n'")
            f.write("end\n")
            f.write("if H_3d_ligands_hybridization_mlct == 1 then\n")
            f.write("    Operators = {H_i, Ssqr, Lsqr, Jsqr, Sk, Lk, Jk, Tk, ldots_3d, N_2p, N_3d, N_L2, 'dZ'}\n")
            f.write(r"    header = 'Analysis of the initial Hamiltonian:\n'")
            f.write("\n")
            f.write(
                r"    header = header .. '===========================================================================================================================================\n'")
            f.write("\n")
            f.write(
                r"    header = header .. 'State         <E>     <S^2>     <L^2>     <J^2>      <Sk>      <Lk>      <Jk>      <Tk>     <l.s>    <N_2p>    <N_3d>    <N_L2>          dZ\n'")
            f.write("\n")
            f.write(
                r"    header = header .. '===========================================================================================================================================\n'")
            f.write("\n")
            f.write(
                r"    footer = '===========================================================================================================================================\n'")
            f.write("\n")
            f.write("end\n")
            f.write("io.write(header)\n")
            f.write("for i, Psi in ipairs(Psis_i) do\n")
            f.write("    io.write(string.format('!*! %5d', i))\n")
            f.write("    for j, Operator in ipairs(Operators) do\n")
            f.write("        if j == 1 then\n")
            f.write("            io.write(string.format('%12.6f', Complex.Re(Psi * Operator * Psi)))\n")
            f.write("        elseif Operator == 'dZ' then\n")
            f.write("            io.write(string.format('%12.2E', dZ[i]))\n")
            f.write("        else\n")
            f.write("            io.write(string.format('%10.4f', Complex.Re(Psi * Operator * Psi)))\n")
            f.write("        end\n")
            f.write("    end\n")
            f.write(r"    io.write('\n')")
            f.write("\n")
            f.write("end\n")
            f.write("io.write(footer)\n")
            f.write("if next(spectra) == nil then\n")
            f.write("    return\n")
            f.write("end\n")
            f.write("E_gs_i = Psis_i[1] * H_i * Psis_i[1]\n")
            f.write("if CalculationRestrictions == nil then\n")
            f.write("    Psis_f = Eigensystem(H_f, FinalRestrictions, 1)\n")
            f.write("else\n")
            f.write(
                "    Psis_f = Eigensystem(H_f, FinalRestrictions, 1, {{'restrictions', CalculationRestrictions}})\n")
            f.write("end\n")
            f.write("Psis_f = {" + "Psis_f}\n")
            f.write("E_gs_f = Psis_f[1] * H_f * Psis_f[1]\n")
            f.write(f"Eedge1 = {Edge}\n")
            f.write("DeltaE = E_gs_f - E_gs_i\n")
            f.write(f"Emin = ({Emin1} - Eedge1 ) + DeltaE\n")
            f.write(f"Emax = ({Emax1} - Eedge1 ) + DeltaE\n")
            f.write(f"NE = {NE1}\n")
            f.write(f"Gamma = {Gamma}\n")
            f.write(f"DenseBorder = {DenseBorder}\n")
            f.write("if CalculationRestrictions == nil then\n")
            f.write(
                "    G_2p_3d = CreateSpectra(H_f, T_2p_3d, Psis_i, {{'Emin', Emin}, {'Emax', Emax}, {'NE', NE}, {'Gamma', Gamma}, {'DenseBorder', DenseBorder}})\n")
            f.write("else\n")
            f.write(
                "    G_2p_3d = CreateSpectra(H_f, T_2p_3d, Psis_i, {{'Emin', Emin}, {'Emax', Emax}, {'NE', NE}, {'Gamma', Gamma}, {'restrictions', CalculationRestrictions}, {'DenseBorder', DenseBorder}})\n")
            f.write("end\n")
            f.write("-- Create a list with the Boltzmann probabilities for a given operator\n")
            f.write("-- and state.\n")
            f.write("dZ_2p_3d = {}\n")
            f.write("for i in ipairs(T_2p_3d) do\n")
            f.write("    for j in ipairs(Psis_i) do\n")
            f.write("        table.insert(dZ_2p_3d, dZ[j])\n")
            f.write("    end\n")
            f.write("end\n")
            f.write("Pcl_2p_3d = 2\n")
            f.write("spectrum = 'Isotropic'\n")
            f.write("if ValueInTable(spectrum, spectra) then\n")
            f.write("    Giso = GetSpectrum(G_2p_3d, T_2p_3d, Psis_i, indices_2p_3d[spectrum], dZ_2p_3d)\n")
            f.write("    Giso = Giso / 3 / Pcl_2p_3d\n")
            f.write("    SaveSpectrum(Giso, 'iso')\n")
            f.write("end\n")
            f.write("spectrum = 'Circular Dichroism'\n")
            f.write("if ValueInTable(spectrum, spectra) then\n")
            f.write("    Gr = GetSpectrum(G_2p_3d, T_2p_3d, Psis_i, indices_2p_3d[spectrum][1], dZ_2p_3d)\n")
            f.write("    Gl = GetSpectrum(G_2p_3d, T_2p_3d, Psis_i, indices_2p_3d[spectrum][2], dZ_2p_3d)\n")
            f.write("    Gr = Gr / Pcl_2p_3d\n")
            f.write("    Gl = Gl / Pcl_2p_3d\n")
            f.write("    SaveSpectrum(Gr, 'r')\n")
            f.write("    SaveSpectrum(Gl, 'l')\n")
            f.write("    SaveSpectrum(Gr - Gl, 'cd')\n")
            f.write("end\n")
            f.write("spectrum = 'Linear Dichroism'\n")
            f.write("if ValueInTable(spectrum, spectra) then\n")
            f.write("    Gv = GetSpectrum(G_2p_3d, T_2p_3d, Psis_i, indices_2p_3d[spectrum][1], dZ_2p_3d)\n")
            f.write("    Gh = GetSpectrum(G_2p_3d, T_2p_3d, Psis_i, indices_2p_3d[spectrum][2], dZ_2p_3d)\n")
            f.write("    Gv = Gv / Pcl_2p_3d\n")
            f.write("    Gh = Gh / Pcl_2p_3d\n")
            f.write("    SaveSpectrum(Gv, 'v')\n")
            f.write("    SaveSpectrum(Gh, 'h')\n")
            f.write("    SaveSpectrum(Gv - Gh, 'ld')\n")
            f.write("end\n")

    def run(self):
        """
        Runs Quanty with the input file specified by Label.lua, and
        returns the standard output and error (if any)

        Arguments:
          Label    : Name of the system (string)
          Qty_path : Path to the Quanty executable (string)

        Returns:
          result   : a subprocess CompltedProcess object
                     result.stdout  is the standard output
                     result.stderr  is the standard error
        """
        command = [self.Qty_path, self.filename]
        result = subprocess.run(command, capture_output=True, text=True)
        self.result = result

    def treat_output(self):
        """
        From the standard output of a Quanty calculation with the XAS_Template, 
        it extracts the relevant expctation value

        Arguments:
          Rawout   : a subprocess CompltedProcess object

        Returns:
          A dictionary with the relevant expectation values

        """
        if not self.result:
            print('Run calculation first!')
            return

        Rawout = self.result
        out = Rawout.stdout.split('\n')
        rline = 0
        for iline in range(len(out)):
            if '!*!' in out[iline]:
                # print(out[iline-2])
                # print(out[iline])
                rline = iline

        Odata = out[rline].split()
        E = Odata[2]
        S2 = Odata[3]
        L2 = Odata[4]
        J2 = Odata[5]
        S_k = Odata[6]
        L_k = Odata[7]
        J_k = Odata[8]
        T_k = Odata[9]
        LdotS = Odata[10]
        self.outdic = {'E': E, 'S2': S2, 'L2': L2, 'J2': J2, 'S_k': S_k, 'L_k': L_k, 'J_k': J_k, 'T_k': T_k,
                       'LdotS': LdotS}

    def post_proc(self):

        import numpy as np
        from scipy.integrate import romb
        from IPython.display import display
        from tabulate import tabulate
        import matplotlib.pyplot as plt

        if not self.outdic:
            self.treat_output()

        outdic = self.outdic
        label = self.ion + '_XAS'
        ion = self.ion
        xz = np.loadtxt(os.path.join(self.path, label + '_iso.spec'), skiprows=5)
        mcd = np.loadtxt(os.path.join(self.path, label + '_cd.spec'), skiprows=5)
        xl = np.loadtxt(os.path.join(self.path, label + '_l.spec'), skiprows=5)
        xr = np.loadtxt(os.path.join(self.path, label + '_r.spec'), skiprows=5)
        mld = np.loadtxt(os.path.join(self.path, label + '_ld.spec'), skiprows=5)
        xv=np.loadtxt(os.path.join(self.path, label + '_v.spec'), skiprows=5)
        xh=np.loadtxt(os.path.join(self.path, label + '_h.spec'), skiprows=5)
        mcd2 = xr.copy()
        mcd2[:, 2] = xl[:, 2] - xr[:, 2]
        npts = np.shape(xz)[0]
        mcd2 = mcd
        # TOTAL spectra
        xas = xz.copy()
        xas[:, 2] = xz[:, 2] + xl[:, 2] + xr[:, 2]

        xas0 = xz.copy()
        xas0[:, 2] = (xl[:, 2] + xr[:, 2]) / 2 + xl[:, 2] + xr[:, 2]

        dx = xz.copy()
        dx[:, 2] = xl[:, 2] + xr[:, 2] - 2 * xz[:, 2]
        # ### Integration using Trapezoidal rule

        nh = 10 - self.params['Nelec']

        trapz = False
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

            lz = -2 * nh * np.trapz(mcd2[:, 2], mcd2[:, 0]) / tot
            szef = 3 / 2 * nh * (
                    np.trapz(mcd2[0:npts // 2, 2], mcd2[0:npts // 2, 0]) -
                    2 * np.trapz(mcd2[npts // 2:, 2], mcd2[npts // 2:, 0])
            ) / tot
            lz0 = -2 * nh * np.trapz(mcd2[:, 2], mcd2[:, 0]) / tot0
            szef0 = 3 / 2 * nh * (
                    np.trapz(mcd2[0:npts // 2, 2], mcd2[0:npts // 2, 0]) -
                    2 * np.trapz(mcd2[npts // 2:, 2], mcd2[npts // 2:, 0])
            ) / tot0

        else:
            print(len(mcd2[npts // 2:, 2]), len(mcd2[0:npts // 2 + 1]))
            mydelta = mcd2[1, 0] - mcd2[0, 0]
            lz = -2 * nh * romb(mcd2[:, 2], float(mydelta)) / tot
            szef = 3 / 2 * nh * (
                    romb(mcd2[0:npts // 2 + 1, 2], float(mydelta)) -
                    2 * romb(mcd2[npts // 2:, 2], float(mydelta))
            ) / tot
            lz0 = -2 * nh * romb(mcd2[:, 2], float(mydelta)) / tot0
            szef0 = 3 / 2 * nh * (
                    romb(mcd2[0:npts // 2 + 1, 2], float(mydelta)) -
                    2 * romb(mcd2[npts // 2:, 2], float(mydelta))
            ) / tot0

        tfmt = 'html'

        Lz_t = float(outdic['L_k'])
        Sz_t = float(outdic['S_k'])
        Tz_t = float(outdic['T_k'])
        Seff_t = float(outdic['S_k']) + float(outdic['T_k'])

        table1 = [[r'L$_z$', r'T$_{z}$',  r'S$_{z}$', r'S$_{eff}$'],
                  [Lz_t, Tz_t, Sz_t, Seff_t]]

        print("Theoretical values (Quanty):")
        display(tabulate(table1, headers='firstrow', tablefmt=tfmt))

        print("Sum rules :")
        table2 = [[r'sL$_z$', 'sS$_{eff}$'], [lz, szef]]
        display(tabulate(table2, headers='firstrow', tablefmt=tfmt))

        print("Sum rules 0:")
        table3 = [[r's$_0$L$_z$', 's$_0$S$_{eff}$'], [lz0, szef0]]
        display(tabulate(table3, headers='firstrow', tablefmt=tfmt))

        print("Deviations:")
        table4 = [[r'$\Delta$XAS (%)', r'$\Delta$L$_{z}$ (%)', r'$\Delta$S$_{eff}$ (%)',
                   r'$\Delta_0$L$_{z}$ (%)', r'$\Delta_0$S$_{eff}$ (%)'],
                  [deltaXas * 100, 100 * (abs(Lz_t) - abs(lz)) / Lz_t,
                   100 * (abs(Seff_t) - abs(szef)) / Seff_t,
                   100 * (abs(Lz_t) - abs(lz0)) / Lz_t,
                   100 * (abs(Seff_t) - abs(szef0)) / Seff_t]]
        display(tabulate(table4, headers='firstrow', tablefmt=tfmt))

        element_data = self.xdat['elements'][self.ion]['charges'][self.charge]
        iondata = element_data['symmetries'][self.symm]['experiments']['XAS']['edges']
        edge = iondata['L2,3 (2p)']['axes'][0][4]

        figs, ax = plt.subplots(3, sharex=True, figsize=(10, 12))
        figs.suptitle(label, fontsize=16)
        ax[0].plot(xz[:, 0] + edge, xz[:, 2], 'r--', label='z-pol')
        ax[0].plot(xl[:, 0] + edge, xl[:, 2], 'b', label='left')
        ax[0].plot(xr[:, 0] + edge, xr[:, 2], 'g', label='right')
        ax[0].set_xlim(-10 + edge, 20 + edge)
        ax[0].legend(fontsize=14)
        ax[0].set_ylabel('Intensity (a.u.)', fontsize=16)

        ax[1].plot(xas[:, 0] + edge, xas[:, 2] / 3, 'k', label='average')
        ax[1].plot(mcd[0:npts // 2, 0] + edge, mcd[0:npts // 2, 2], 'r', label=r'L$_3$')
        ax[1].plot(mcd[npts // 2:, 0] + edge, mcd[npts // 2:, 2], 'b', label=r'L$_2$')
        ax[1].set_xlim(-10 + edge, 20 + edge)
        ax[1].legend(fontsize=14, title='XMCD')
        ax[1].set_ylabel('Intensity (a.u.)', fontsize=16)
        ax[1].set_xlabel('Energy (eV)', fontsize=16)

        ax[2].plot(xas[:, 0] + edge, xas[:, 2] / 3, 'k', label='average')
        ax[2].plot(xv[: , 0] + edge, xv[:, 2], 'b--', label=r'Vertical polarization')
        ax[2].plot(xh[: , 0] + edge, xh[:, 2], 'g:', label=r'Horizontal polarization')
        ax[2].plot(mld[: , 0] + edge, mld[:, 2], 'r', label=r'Linear dichroism1')
        ax[2].plot(mld[: , 0] + edge, xv[:,2]-xh[:,2], 'k:', label=r'Linear dichroism2')
        ax[2].set_xlim(-10 + edge, 20 + edge)
        ax[2].legend(fontsize=14, title='XMLD')
        ax[2].set_ylabel('Intensity (a.u.)', fontsize=16)
        ax[2].set_xlabel('Energy (eV)', fontsize=16)

        plt.show()

    def post_proc_output_only(self):

        import numpy as np

        label = self.ion + '_XAS'
        xz = np.loadtxt(os.path.join(self.path, label + '_iso.spec'), skiprows=5)
        mcd = np.loadtxt(os.path.join(self.path, label + '_cd.spec'), skiprows=5)
        xl = np.loadtxt(os.path.join(self.path, label + '_l.spec'), skiprows=5)
        xr = np.loadtxt(os.path.join(self.path, label + '_r.spec'), skiprows=5)

        mcd2 = xr.copy()
        mcd2[:, 2] = xl[:, 2] - xr[:, 2]

        # TOTAL spectra
        xas = xz.copy()
        xas[:, 2] = xz[:, 2] + xl[:, 2] + xr[:, 2]

        xas0 = xz.copy()
        xas0[:, 2] = (xl[:, 2] + xr[:, 2]) / 2 + xl[:, 2] + xr[:, 2]

        dx = xz.copy()
        dx[:, 2] = xl[:, 2] + xr[:, 2] - 2 * xz[:, 2]

        element_data = self.xdat['elements'][self.ion]['charges'][self.charge]
        iondata = element_data['symmetries'][self.symm]['experiments']['XAS']['edges']
        edge = iondata['L2,3 (2p)']['axes'][0][4]

        output = {
            'zpol_energy': xz[:, 0] + edge,
            'zpol_xas': xz[:, 2],
            'xas_left_energy': xl[:, 0] + edge,
            'xas_left': xl[:, 2],
            'xas_right_energy': xr[:, 0] + edge,
            'xas_right': xr[:, 2],
            'average_energy': xas[:, 0] + edge,
            'average': xas[:, 2] / 3,
            'xmcd_energy': mcd[:, 0] + edge,
            'xmcd': mcd[:, 2],
        }
        return output

    def run_all_with_output(self):
        self.write_header()
        self.H_init()
        self.setH_terms()
        self.set_electrons()
        self.define_atomic_term()
        self.define_crystal_field_term()
        self.define_external_field_term()
        self.setTemperature()
        self.setRestrictions()
        self.set_iterative_solver()
        self.set_spectra_functions()
        self.define_transitions([0, 0, 1], [0, 1, 0], [1, 0, 0])
        self.set_spectra_lists()
        self.calculate_and_save_spectra()
        self.run()  # run quanty!
        return self.post_proc_output_only()
    
