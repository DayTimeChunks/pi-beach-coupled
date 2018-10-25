import numpy as np
import csv
# from scipy.stats.distributions import norm
#
# from pcraster._pcraster import *
# from pcraster.framework import *
# from model_v2 import *

# Option to use pyDOE (same sampling results)
# https://pythonhosted.org/pyDOE/randomized.html#randomized
# from pyDOE import *

print("Using bounds 16!")

def scale(bounds):
    """
    :param bounds: parameter bounds to test
    :return: normalized bounds by the upper bound,
        that means that all maximum values are == 1,
        so the bounds need to be re-scaled inside BEACH.
    """
    scaled = []
    for e in bounds:
        se = []
        up = max(e[0], e[1])
        se.append(e[0] / up)
        se.append(e[1] / up)
        scaled.append(se)

    return scaled


def get_upper(bounds):
    upper = []
    for e in bounds:
        upper.append(e[1])
    return upper


def saveLHSmatrix(param_values):
    np.savetxt("tss\\output\\lhs_vectors.txt", param_values)


def get_runs(params):
    runs = int(params.shape[0])
    print("Runs: ", int(runs))
    return runs


def get_problem():
    # Will be scaled with scale()
    bounds = [[0.85, 0.99],  # z3_factor
              [0.01, 1.],  # 'cZ0Z1'
              [0.2, 0.6],  # 'cZ'
              [0.01, 1.],  # cadr
              [1500.0, 3650.0],  # k_g
              [0.01, 1.],  # gamma01,
              [0.01, 1.],  # gammaZ
              [0.1, 1.],  # f_transp
              [0.1, 0.9],  # f_evap
              [0.01, 0.05],  # f_oc,
              [0.3, 1700],  # k_oc, -> New: max Kd = 85
              [0.01, 0.5],  # beta_runoff
              [150.0, 3000.0],  # age_rate, -> New: min was 10.0 (EST)
              [85.0, 350.0],  # dt_50_ab,  -> New: min was 65.0 (EST)
              [13.0, 23.0],  # dt_50_ref -> New: EST
              [1.7, 3.2],  # epsilon (in absolute, convert to negative!!) -> New: EST paper
              [0.01, 1.0]]  # beta_moisture

    names = ['z3_factor',
             'cZ0Z1', 'cZ',
             'c_adr',
             'k_g',
             'gamma01', 'gammaZ',
             'f_transp',
             'f_evap',
             'f_oc', 'k_oc',
             'beta_runoff',
             'dt_50_aged',
             'dt_50_ab',
             'dt_50_ref',
             'epsilon_iso',
             'beta_moisture'
             ]

    problem = {
        'num_vars': len(bounds),
        'names': names,
        'bounds': scale(bounds),
        'upper': get_upper(bounds)
    }
    return problem


def get_vector_test():
    names = ['z3_factor',
             'cZ0Z1', 'cZ',
             'c_adr',
             'k_g',
             'gamma01', 'gammaZ',
             'f_transp',
             'f_evap',
             'f_oc', 'k_oc',
             'beta_runoff',
             'dt_50_aged',
             'dt_50_ab',
             'dt_50_ref',
             'epsilon_iso',
             'beta_moisture'
             ]
    # print("names length : " + str(len(names)))

    ini_path = 'csv\\initial.csv'
    test_param = dict()  # Dictionary to store the values
    with open(ini_path, 'r') as f:
        reader = csv.reader(f, delimiter=',')
        for row in reader:
            test_param[row[0].strip()] = float(row[1])
    name_values = []
    for i in range(len(names)):
        name_values.append(test_param[names[i]])
    # print("vales length : " + str(len(name_values)))
    return name_values
