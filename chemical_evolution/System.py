import numpy as np
from scipy.integrate import ode


class ChemEvo:
    def __init__(self, all_species, all_reactions, temp,
                 init_t=0, integrator='vode', **kwargs):
        self.all_species = all_species
        self.all_reactions = all_reactions
        self.rhs_dict = {species: species.fetch_rhs() for species in self.all_species}

        for rxn in self.all_reactions:
            rxn.calc_rate(temp)
        self.t = init_t

        def sys_rhs(t, y):
            state = {**{rxn.rate_symbol: rxn.rate_val for rxn in self.all_reactions},
                     **{species.symbol: species.amount for species in self.all_species}}
            deltas = [
                self.rhs_dict[species].subs(**state)
                for species in self.all_species
            ]
            return deltas

        init_y = [species.amount for species in all_species]
        self.integrator = ode(sys_rhs).set_integrator(integrator, **kwargs)
        self.integrator.set_initial_value(init_y, self.t)

    def update(self, dt):
        updated_amounts = self.integrator.integrate(self.integrator.t + dt)
        for amt, species in zip(updated_amounts, self.all_species):
            species.amount_hist.append(amt)
        self.t += dt
