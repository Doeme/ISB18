
import cobra.test
import os
from os.path import join
import itertools
import pickle
import sys
import math




filenames = [
    'results_question_1_acetate.dat',
    'results_question_1_d-lactate.dat',
    'results_question_1_ethanol.dat',
    'results_question_1_l-lactate.dat',
    'results_question_1_mutial_lactate.dat',
    'results_question_1_succinate.dat']

for filename in filenames:
    with open(filename, 'rb') as pfile:
        name = pickle.load(pfile)
        #data = []
        #while True:
        #    try:
        #        data.append(pickle.load(pfile))
        #    except (EOFError, pickle.UnpicklingError):
        #        break
        data = pickle.load(pfile)

        status = [ item[1] for item in data ]
        cnt = status.count('threshold')
        cnt = len(data)

        print("{}\t{}".format(filename, cnt))



print("END")
