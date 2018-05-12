import numpy as np
import sympy as sm


class Species:
    """
    This class models each species
    """
    def __init__(self, name, init_amount, symbol=None):
        self.name = name
        if symbol is not None:
            self.symbol = symbol
        else:
            self.symbol = sm.Symbol(name)
        self.amount = init_amount
        self.amount_hist = [init_amount]

    def _find_formation(self, all_reactions):
        """
        From the list of reactions passed as an argument, find those that are
        constructive for this species.
        """
        return [rxn for rxn in all_reactions if self.symbol in rxn.rhs.atoms()]

    def _find_destruction(self, all_reactions):
        """
        From the list of reactions passed as an argument, find those that are
        destructive for this species.
        """
        return [rxn for rxn in all_reactions if self.symbol in rxn.lhs.atoms()]

    def fetch_rhs(self, all_reactions):
        """
        Solve the right hand equation for the $dS/dt$ of this species, from the
        list of reactions passed as an argument.
        """
        dSdt = 0
        for rxn in self._find_formation(all_reactions):
            term = rxn.rate_symbol
            for atom in rxn.lhs.atoms():
                term *= atom
            dSdt += term
        for rxn in self._find_destruction(all_reactions):
            term = - rxn.rate_symbol
            for atom in rxn.lhs.atoms():
                term *= atom
            dSdt += term
        return dSdt

