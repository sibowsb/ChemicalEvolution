import sympy as sm
import numpy as np


class Reaction:
    def __init__(self, lhs, rhs, rate=0, rate_symbol=None):
        self.lhs = lhs
        self.rhs = rhs
        if rate_symbol is None:
            self.rate_symbol = sm.Symbol('k')
        else:
            self.rate_symbol = rate_symbol
        self.rate_val = rate
