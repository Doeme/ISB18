
import cobra.test
import os
from os.path import join
import itertools
import pickle
import sys
import math






# Load data
with open('results_question_1_EX_ac(e).dat', 'rb') as f:
    ex_reac = pickle.load(f)
    data = pickle.load(f)

for elem in data:
    print("{}\t{}\t{}".format(elem[1], elem[2], elem[3]))



print("END")
