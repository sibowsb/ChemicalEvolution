import numpy as np
import sympy as sm
from re import findall


def calc_rate(T, rx):
    T_eV = T / 11605.
    log_T = np.log(T)
    log_T_eV = np.log(T_eV)
    if rx == 0:
        rate = np.exp(- 32.71396786375
                      + 13.53655609057 * log_T_eV
                      - 5.739328757388 * log_T_eV ** 2
                      + 1.563154982022 * log_T_eV ** 3
                      - 0.2877056004391 * log_T_eV ** 4
                      + 0.03482559773736999 * log_T_eV ** 5
                      - 0.00263197617559 * log_T_eV ** 6
                      + 0.0001119543953861 * log_T_eV ** 7
                      - 2.039149852002e-6 * log_T_eV ** 8)
    elif rx == 1:
        rate = 4.881357e-6 * T ** (-1.5) * (1. + 1.14813e2 * T ** (-0.407)) ** (-2.242)
    elif rx == 2:
        rate = np.exp(- 44.09864886561001
                      + 23.91596563469 * log_T_eV
                      - 10.75323019821 * log_T_eV ** 2
                      + 3.058038757198 * log_T_eV ** 3
                      - 0.5685118909884001 * log_T_eV ** 4
                      + 0.06795391233790001 * log_T_eV ** 5
                      - 0.005009056101857001 * log_T_eV ** 6
                      + 0.0002067236157507 * log_T_eV ** 7
                      - 3.649161410833e-6 * log_T_eV ** 8)
    elif rx == 3:
        rate = (1.54e-9 * (1. + 0.3 / np.exp(8.099328789667 / T_eV)) /
                (np.exp(40.49664394833662 / T_eV) * T_eV ** 1.5)
                + 3.92e-13 / T_eV ** 0.6353)
    elif rx == 4:
        rate = np.exp(-68.71040990212001
                      + 43.93347632635 * log_T_eV
                      - 18.48066993568 * log_T_eV ** 2
                      + 4.701626486759002 * log_T_eV ** 3
                      - 0.7692466334492 * log_T_eV ** 4
                      + 0.08113042097303 * log_T_eV ** 5
                      - 0.005324020628287001 * log_T_eV ** 6
                      + 0.0001975705312221 * log_T_eV ** 7
                      - 3.165581065665e-6 * log_T_eV ** 8)

    elif rx == 5:
        rate = 7.8155e-5 * T ** (-1.5) * (1. + 2.0189e2 * T ** (-0.407)) ** (-2.242)

    elif rx == 6:
        rate = 3.0e-16 * (T / 3e2) ** 0.95 * np.exp(-T / 9.32e3)

    elif rx == 7:
        rate = 1.35e-9 * (T ** 9.8493e-2 + 3.2852e-1 * T ** 5.5610e-1 + 2.771e-7 * T ** 2.1826) \
               / (1. + 6.191e-3 * T ** 1.0461 + 8.9712e-11 * T ** 3.0424 + 3.2576e-14 * T ** 3.7741)

    elif rx == 8:
        rate = 2.10e-20 * (T / 30.0) ** (-0.15)

    elif rx == 9:
        rate = 6.0e-10 * np.ones_like(T)

    elif rx == 10:
        rate = (np.exp(-21237.15 / T)
                * (- 3.3232183e-07
                   + 3.3735382e-07 * log_T
                   - 1.4491368e-07 * log_T ** 2
                   + 3.4172805e-08 * log_T ** 3
                   - 4.7813720e-09 * log_T ** 4
                   + 3.9731542e-10 * log_T ** 5
                   - 1.8171411e-11 * log_T ** 6
                   + 3.5311932e-13 * log_T ** 7))

    elif rx == 11:
        rate = 4.4886e-9 * T ** 0.109127 * np.exp(-101858. / T)

    elif rx == 12:
        rate = 1.0670825e-10 * T_eV ** 2.012 / (np.exp(4.463 / T_eV) * (1. + 0.2472 * T_eV) ** 3.512)

    elif rx == 13:
        rate = np.exp(-18.01849334273
                      + 2.360852208681 * log_T_eV
                      - 0.2827443061704 * log_T_eV ** 2
                      + 0.01623316639567 * log_T_eV ** 3
                      - 0.03365012031362999 * log_T_eV ** 4
                      + 0.01178329782711 * log_T_eV ** 5
                      - 0.001656194699504 * log_T_eV ** 6
                      + 0.0001068275202678 * log_T_eV ** 7
                      - 2.631285809207e-6 * log_T_eV ** 8)

    elif rx == 14:
        rate = np.exp(-20.37260896533324
                      + 1.139449335841631 * log_T_eV
                      - 0.1421013521554148 * log_T_eV ** 2
                      + 0.00846445538663 * log_T_eV ** 3
                      - 0.0014327641212992 * log_T_eV ** 4
                      + 0.0002012250284791 * log_T_eV ** 5
                      + 0.0000866396324309 * log_T_eV ** 6
                      - 0.00002585009680264 * log_T_eV ** 7
                      + 2.4555011970392e-6 * log_T_eV ** 8
                      - 8.06838246118e-8 * log_T_eV ** 9)

    elif rx == 15:
        rate = 2.4e-6 * (1. + T / 2e4) / np.sqrt(T)

    elif rx == 16:
        rate = 1e-8 * T ** (-0.4)

    elif rx == 17:
        rate = 1e-8 * np.ones_like(T)

    elif rx == 18:
        rate = 5e-7 * np.sqrt(100. / T)

    else:
        rate = 0 * np.ones_like(T)

    return rate


def parse_reactions():
    source = r'''- ${\rm H + e^{-}} \rightarrow {\rm H^{+} + e^{-} + e^{-}}$
                 - ${\rm H^{+} + e^{-}} \rightarrow {\rm H + \gamma}$
                 - ${\rm He + e^{-}} \rightarrow {\rm He^{+} + e^{-} + e^{-}}$
                 - ${\rm He^{+} + e^{-}} \rightarrow {\rm He + \gamma}$
                 - ${\rm He^{+} + e^{-}} \rightarrow {\rm He^{++} + e^{-} + e^{-}}$
                 - ${\rm He^{++} + e^{-}} \rightarrow {\rm He^{+} + \gamma}$
                 - ${\rm H + H} \rightarrow {\rm H^{+} + e^{-} + H}$
                 - ${\rm H + He} \rightarrow {\rm H^{+} + e^{-} + He}$
                 - ${\rm H + \gamma} \rightarrow {\rm H^{+} + e^{-}}$
                 - ${\rm He + \gamma} \rightarrow {\rm He^{+} + e^{-}}$
                 - ${\rm He^{+} + \gamma} \rightarrow {\rm He^{++} + e^{-}}$
                 - ${\rm H + e^{-}} \rightarrow {\rm H^{-} + \gamma}$
                 - ${\rm H^{-} + H} \rightarrow {\rm H_{2} + e^{-}}$
                 - ${\rm H + H^{+}} \rightarrow {\rm H_{2}^{+} + \gamma}$
                 - ${\rm H_{2}^{+} + H} \rightarrow {\rm H_{2} + H^{+}}$
                 - ${\rm H_{2} + H^{+}} \rightarrow {\rm H_{2}^{+} + H}$
                 - ${\rm H_{2} + e^{-}} \rightarrow {\rm H + H + e^{-}}$
                 - ${\rm H_{2} + H} \rightarrow {\rm H + H + H}$
                 - ${\rm H^{-} + e^{-}} \rightarrow {\rm H + e^{-} + e^{-}}$
                 - ${\rm H^{-} + H} \rightarrow {\rm H + e^{-} + H}$
                 - ${\rm H^{-} + H^{+}} \rightarrow {\rm H + H}$
                 - ${\rm H^{-} + H^{+}} \rightarrow {\rm H_{2}^{+} + e^{-}}$
                 - ${\rm H_{2}^{+} + e^{-}} \rightarrow {\rm H + H}$
                 - ${\rm H_{2}^{+} + H^{-}} \rightarrow {\rm H_{2} + H}$
                 - ${\rm H + H + H} \rightarrow {\rm H_{2} + H}$
                 - ${\rm H + H + H_{2}} \rightarrow {\rm H_{2} + H_{2}}$
                 - ${\rm H^{-} + \gamma} \rightarrow {\rm H + e^{-}}$
                 - ${\rm H_{2}^{+} + \gamma} \rightarrow {\rm H + H^{+}}$
                 - ${\rm H_{2} + \gamma} \rightarrow {\rm H_{2}^{+} + e^{-}}$
                 - ${\rm H_{2}^{+} + \gamma} \rightarrow {\rm H^{+} + H^{+} + e^{-}}$
                 - ${\rm H_{2} + \gamma} \rightarrow {\rm H + H}$
                 - ${\rm H + H + grain} \rightarrow {\rm H_{2} + grain}$'''

    formulas = findall(r'\${\\rm (.*)} \\rightarrow {\\rm (.*)}\$',
                       source.replace('H}', 'H }').replace('He}', 'He }'))

    replace_rules = [
        ('H^{-}', 'HM'), ('H^{+}', 'HII'),
        ('H_{2}^{+}', 'H2II'), ('H_{2}', 'H2I'), ('H ', 'HI '),
        ('He^{++}', 'HeIII'), ('He^{+}', 'HeII'), ('He ', 'HeI '),
        ('e^{-}', 'de'), ('\\gamma', 'gma'), ('grain', 'grn')
    ]

    formulas_clean = []
    for i in range(len(formulas)):
        lhs, rhs = formulas[i]
        for latex, plain in replace_rules:
            lhs = lhs.replace(latex, plain)
            rhs = rhs.replace(latex, plain)
        formulas_clean.append((lhs, rhs))

    return formulas_clean

    # ks = sm.sympify(
    #     'k1, k2, k3, k4, k5, k6, k7, k8, k9, k10, k11, k12, '
    #     'k13, k14, k15, k16, k17, k18, k19, k20, k21, k22, k23, k24, '
    #     'k25, k26, k27, k28, k29, k30, k31, k32'
    # )
    #
    # all_reactions = [
    #     eval('({lhs}), ({rhs}), ks[{i}]'.format(lhs=lhs, rhs=rhs, i=i))
    #     for i, (lhs, rhs) in enumerate(formulas_clean)
    # ]
    #
    # return all_reactions
