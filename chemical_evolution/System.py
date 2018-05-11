import numpy as np
from scipy.integrate import ode


class ChemEvo:
    def __init__(self, all_species, all_reactions,
                 init_t=0, integrator='vode', **kwargs):
        self.all_species = all_species
        self.all_reactions = all_reactions
        self.rhs_dict = {species: species.fetch_rhs(all_reactions) for species in self.all_species}

        self.t = init_t
        
        def sys_rhs(t, y):
            states = {**{rxn.rate_symbol: rxn.rate_val for rxn in self.all_reactions},
                      **{species.symbol: y[i] for i, species in enumerate(self.all_species)}}
            delta = [self.rhs_dict[species].subs(states) for species in self.all_species]
            return delta

        self.sys_rhs = sys_rhs

        init_y = [species.amount for species in all_species]
        self.y = init_y
        self.integrator = ode(self.sys_rhs)
        self.integrator.set_integrator(integrator, **kwargs)
        self.integrator.set_initial_value(init_y, self.t)

    def update(self, dt):
        updated_amounts = self.integrator.integrate(self.integrator.t + dt)
        for amt, species in zip(updated_amounts, self.all_species):
            species.amount_hist.append(amt)
        self.y = updated_amounts
        self.t += dt

        
def sys_rhs(y, t, sys):
    states = {**{rxn.rate_symbol: rxn.rate_val for rxn in sys.all_reactions},
              **{species.symbol: y[i] for i, species in enumerate(sys.all_species)}}
    delta = [sys.rhs_dict[species].subs(states) for species in sys.all_species]
    return delta