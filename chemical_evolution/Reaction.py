import sympy as sm
import numpy as np


class Reaction:
    def __init__(self, lhs, rhs):
        self.lhs = lhs
        self.rhs = rhs
        self.rate_symbol = sm.Symbol('k')
        self.rate_val = np.nan

    def calc_rate(self, temp):
        self.rate_val = 0
        return self.rate_val