import numpy as np
import itertools
import os

def test_reference_cases(dirname, test_fun):
    """
    Run several tests by combining PCV parameters and compare outcomes with
    reference data from R.
    """

    tests_path = '../.tests/'
    data_path = tests_path + 'data/'
    refval_path = tests_path + dirname + '/'

    # check if reference files exist
    if  not os.path.exists(tests_path) or \
        not os.path.exists(data_path) or \
        not os.path.exists(data_path + 'corn.csv') or \
        not os.path.exists(refval_path):

        print('can not find reference files, skipping.')
        return

    # red the dataset for tests
    D = np.genfromtxt('../.tests/data/corn.csv', delimiter=',')
    X = D[:, 1:]
    Y = D[:, :1]

    cv_cases = [{'type':'loo'}, {'type':'ven', 'nseg': 4}, {'type':'ven', 'nseg': 10}]
    ncomp_cases = [1, 10, 20, 30]
    scale_cases = [True, False]

    # loop over all combinations of the parameters
    all_cases = list(itertools.product(ncomp_cases, cv_cases, scale_cases))
    for ncomp, cv, scale in all_cases:

        # file suffix for reference values - depends on all parameters
        file_suffix = '-' + str(ncomp) + '-' + str(scale).upper() + '-' + \
            (cv['type'] + str(cv['nseg']) if 'nseg' in cv else cv['type']) + '.csv'

        # run test
        test_fun(X, Y, ncomp, cv = cv, center = True, scale = scale, path = refval_path, file_suffix = file_suffix)