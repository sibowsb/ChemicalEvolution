#!/usr/bin/env python3

import argparse
import chemical_evolution as chemevo
import nine_species_network as net9
import numpy as np
import matplotlib.pyplot as plt

"""
Initial H ionization fraction
Initial He ionization fraction
Molecular hydrogen fraction
Density
Time to evolve to
Integrator to use inside scipy.integrate.ode

"""

parser = argparse.ArgumentParser(
    description='9-species network simulator. ' \
                'Example: python 9_species_net.py 1 1 1 1 100000 1000000 vode'
)
parser.add_argument('h_frac', help='initial H ionization fraction', type=float)
parser.add_argument('he_frac', help='initial He ionization fraction', type=float)
parser.add_argument('h2_frac', help='initial molecular hydrogen fraction', type=float)
parser.add_argument('density', help='density', type=float)
parser.add_argument('temperature', help='temperature', type=float)
parser.add_argument('end_time', help='time to evolve to', type=float)
parser.add_argument('integrator', help='integrator to use inside scipy.integrate.ode', type=str)
args = parser.parse_args()

# Calculate initial condition
h_amt = args.density * args.h_frac
he_amt = args.density * args.he_frac
h2_amt = args.density * args.h2_frac

# Load species
species_info = net9.LoadingHook.load_species()
all_species = [
    chemevo.Species(name, 1, symbol=symbol)
    for name, symbol in species_info
]
all_species[1].amount = h_amt
all_species[2].amount = he_amt
all_species[4].amount = h2_amt
HI, HII, HM, HeI, HeII, HeIII, H2I, H2II, de, gma, grn = [species.symbol for species in all_species]

# Load reactions
reaction_info = net9.LoadingHook.load_reactions()
all_reactions = [
    chemevo.Reaction(eval(lhs), eval(rhs), rate=net9.HardCoded.calc_rate(1e5, i), rate_symbol=rate_symbol)
    for i, ((lhs, rhs), rate_symbol) in enumerate(reaction_info)
]

# Construct the system
system = chemevo.ChemEvo(all_species, all_reactions)

# Run the simulation
y_hist = [system.y]
for i in range(10000):
    # if i % int(args.end_time / 10) == 0:
    #     print('Progress: %d%%' % (100 * system.t / args.end_time))
    system.update(args.end_time / 10000)
    y_hist.append(system.y)

# Plot
t = np.linspace(0, args.end_time, 10001)
for i in range(11):
    plt.plot(t, [x[i] for x in y_hist], label=system.all_species[i].name)
plt.legend([species.name for species in system.all_species])
plt.xlabel('Time (s)')
plt.ylabel('Fraction')

plt.show()