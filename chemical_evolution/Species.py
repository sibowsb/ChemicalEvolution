import numpy as np
import sympy as sm


class Species:
    def __init__(self, name, init_amount):
        self.symbol = sm.Symbol(name)
        self.amount = init_amount
        self.amount_hist = [init_amount]

    def _find_formation(self, all_reactions):
        return [rxn for rxn in all_reactions if self.symbol in rxn.rhs.atoms()]

    def _find_destruction(self, all_reactions):
        return [rxn for rxn in all_reactions if self.symbol in rxn.lhs.atoms()]

    def fetch_rhs(self, all_reactions):
        dSdt = 0
        for rxn in self._find_formation(all_reactions):
            term = rxn.rate_symbol
            for atom in rxn.lhs.atoms():
                term *= atom.symbol
            dSdt += term
        for rxn in self._find_destruction(all_reactions):
            term = - rxn.rate_symbol
            for atom in rxn.lhs.atoms():
                term *= atom.symbol
            dSdt += term
        return dSdt

