import sympy as sm
from nine_species_network import *


def load_species():
    """
    Load the species info as an iterable of name-symbol pairs, where name is the name of the
    species as a string, and symbol is the sympy symbol for its fraction.
    """
    all_species = sm.sympify(
        'HI, HII, HM, HeI, HeII, HeIII, H2I, H2II, de, gma, grn'
    )
    all_names = ('HI', 'HII', 'HM', 'HeI', 'HeII', 'HeIII', 'H2I', 'H2II', 'de', 'gma', 'grn')
    return zip(all_names, all_species)


def load_reactions():
    """
    Load the reaction info as an iterable of equation-rate pairs, where equation is the LHS 
    and the RHS as a tuple, and rate is the sympy symbol for the reaction rate.
    """
    all_reactions = HardCoded.parse_reactions()
    ks = sm.sympify(
        'k0, k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, '
        'k13, k14, k15, k16, k17, k18, k19, k20, k21, k22, k23, k24, '
        'k25, k26, k27, k28, k29, k30, k31'
    )
    return zip(all_reactions, ks)

