
import cobra.test
import os
from os.path import join
import itertools
import pickle
import sys
import math
import matplotlib.pyplot as plt



filenames = [
#    'results_question_1_acetate.dat',
    'results_question_1_d-lactate.dat',
#    'results_question_1_ethanol.dat',
    'results_question_1_l-lactate.dat',
    'results_question_1_mutial_lactate.dat',
    'results_question_1_succinate.dat']

def load_pickle_file(filename, ftype='lst'):
    with open(filename, 'rb') as pfile:
        name = ''
        data = []
        name = pickle.load(pfile)
        if ftype == 'flat':
            while True:
                try:
                    data.append(pickle.load(pfile))
                except (EOFError, pickle.UnpicklingError):
                    break
        elif ftype == 'lst':
            data = pickle.load(pfile)
        else:
            print("Error: File type \"{}\" is not recognized".format(ftype))

    return (name, data)

def show_all_files():
    for filename in filenames:
        name, data = load_pickle_file(filename)
        
        status = [ item[1] for item in data ]
        cnt = status.count('threshold')
        cnt = len(data)
    
        print("{}\t{}".format(filename, cnt))

def inspect_file(filename, row_delta=20):
    name, data = load_pickle_file(filename)

    proceed = True
    offset = 0
    while proceed:
        start_idx = offset
        end_idx = min(start_idx + row_delta - 1, len(data))
        for idx in range(start_idx, end_idx + 1):
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}".format(idx, data[idx][1], data[idx][2], data[idx][3], data[idx][4], data[idx][5], data[idx][6]))

        if end_idx == len(data):
            print("End of file")
        else:
            user_input = input('Load next {} rows? (y/n):'.format(row_delta))
            yanswers = ['y', 'yes', 'Y', 'Yes', 'YES', 'yup', 'Yup', 'jup', 'Jup',
                    'yep', 'Yep', 'jep', 'Jep', 'mhm', '']
            proceed = user_input in yanswers
            offset = start_idx + row_delta

def draw_histogram_of_file(filename, value='bio'):
    name, data = load_pickle_file(filename)

    val_names = {'bio': 2, 'max': 4, 'min': 6}
    values = [ item[val_names[value]] for item in data ]
    plt.hist(values)
    plt.title('Histogram of column \"{}\" in file {}'.format(value, filename))
    plt.xlabel('values')
    plt.ylabel('frequency')
    plt.show()
        
inspect_file(filenames[0] + '.test')

#draw_histogram_of_file(filenames[0] + '.test', value='bio')
#draw_histogram_of_file(filenames[0] + '.test', value='max')
#draw_histogram_of_file(filenames[0] + '.test', value='min')

print("END")
