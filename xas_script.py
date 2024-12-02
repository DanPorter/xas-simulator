"""
XAS Simulator

Run quanty calculation and perform post-processing for a particular element
"""

from xas_simulator.simulation import create_simulation


# Create and run Quanty Simulation
simulation = create_simulation(
    ion='Ni',  # Ion Name
    ch=2,  # Charge
    beta=0.8,  # Beta
    dq=1.0,  # 10Dq
    mag_field=(0, 0, 0),
    exchange_field=(0, 0, 0.1),
    temperature=1.0,  # T (K)
)
output = simulation.run_all_with_output()
simulation.post_proc()  # create tables and plot

print('Finished!')
