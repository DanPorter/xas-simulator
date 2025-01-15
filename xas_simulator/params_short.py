"""
Atomic & X-Ray parameters
"""

import os
import json

PARAMETERS_FILE = os.path.dirname(__file__) + '/parameters.json'

# x-ray parameters
with open(PARAMETERS_FILE, 'r', encoding='utf-8') as f:
    xray_data = json.load(f, )

# Atomic parameters
parameters = {
    "Cu": {
        "Nelec": 9,
        "zeta_3d": 0.102,
        "F2dd": 12.854,
        "F4dd": 7.980,
        "zeta_2p": 13.498,
        "F2pd": 8.177,
        "G1pd": 6.169,
        "G3pd": 3.510,
        "Xzeta_3d": 0.124,
        "XF2dd": 13.611,
        "XF4dd": 8.457
    },
    "Ni": {
        "Nelec": 8,
        "zeta_3d": 0.083,
        "F2dd": 12.233,
        "F4dd": 7.597,
        "zeta_2p": 11.507,
        "F2pd": 7.720,
        "G1pd": 5.783,
        "G3pd": 3.290,
        "Xzeta_3d": 0.102,
        "XF2dd": 13.005,
        "XF4dd": 8.084
    },
    "Co": {
        "Nelec": 7,
        "zeta_3d": 0.066,
        "F2dd": 11.604,
        "F4dd": 7.209,
        "zeta_2p": 9.748,
        "F2pd": 7.259,
        "G1pd": 5.394,
        "G3pd": 3.068,
        "Xzeta_3d": 0.083,
        "XF2dd": 12.395,
        "XF4dd": 7.707
    },
    "Fe": {
        "Nelec": 6,
        "zeta_3d": 0.052,
        "F2dd": 10.965,
        "F4dd": 6.815,
        "zeta_2p": 8.200,
        "F2pd": 6.792,
        "G1pd": 5.000,
        "G3pd": 2.843,
        "Xzeta_3d": 0.067,
        "XF2dd": 11.778,
        "XF4dd": 7.327
    },
    "Mn": {
        "Nelec": 5,
        "zeta_3d": 0.040,
        "F2dd": 10.315,
        "F4dd": 6.413,
        "zeta_2p": 6.846,
        "F2pd": 6.320,
        "G1pd": 4.603,
        "G3pd": 2.617,
        "Xzeta_3d": 0.053,
        "XF2dd": 11.154,
        "XF4dd": 6.942
    },
     "Cr": {
        "Nelec": 4,
        "zeta_3d": 0.030,
        "F2dd": 9.648,
        "F4dd": 6.001,
        "zeta_2p": 5.668,
        "F2pd": 5.840,
        "G1pd": 4.201,
        "G3pd": 2.387,
        "Xzeta_3d": 0.041,
        "XF2dd": 10.521,
        "XF4dd": 6.551
    },
    "V": {
        "Nelec": 3,
        "zeta_3d": 0.022,
        "F2dd": 8.961,
        "F4dd": 5.576,
        "zeta_2p": 4.650,
        "F2pd": 5.351,
        "G1pd": 3.792,
        "G3pd": 2.154,
        "Xzeta_3d": 0.031,
        "XF2dd": 9.875,
        "XF4dd": 6.152
    },
    "Ti": {
        "Nelec": 2,
        "zeta_3d": 0.016,
        "F2dd": 8.243,
        "F4dd": 5.132,
        "zeta_2p": 3.776,
        "F2pd": 4.849,
        "G1pd": 3.376,
        "G3pd": 1.917,
        "Xzeta_3d": 0.023,
        "XF2dd": 9.21,
        "XF4dd": 5.744
    },
    "Sc": {
        "Nelec": 1,
        "zeta_3d": 0.010,
        "F2dd": 0,
        "F4dd": 0,
        "zeta_2p": 3.032,
        "F2pd": 4.332,
        "G1pd": 2.950,
        "G3pd": 1.674,
        "Xzeta_3d": 0.017,
        "XF2dd": 8.530,
        "XF4dd": 5.321
    }
}