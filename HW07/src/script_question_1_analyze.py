
import cobra.test
import os
from os.path import join
import itertools
import pickle
import sys
import math
import matplotlib.pyplot as plt
from q1_draw_graph import draw_maxmin_by_growth
import numpy as np
import csv


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
    inspect_data(data, fow_delta=row_delta)

def inspect_data(data, row_delta=20):

    proceed = True
    offset = 0
    while proceed and offset < len(data):
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

def draw_histogram_of_file(filename, value='all'):
    name, data = load_pickle_file(filename)
    draw_histogram_of_data(data, value=value)

def draw_histogram_of_data(data, value='all'):
    # Remove all entries where the state was not 'optimal'
    data = [ row for row in data if row[1] == row[3] == row[5] == 'optimal']
    val_names = {'bio': 2, 'max': 4, 'min': 6}
    if value == 'all':
        plt.subplot(2,2,1)
        plt.hist([ item[val_names['bio']] for item in data ])
        plt.title('Histogram of column "bio"')
        plt.xlabel('values')
        plt.ylabel('frequency')

        plt.subplot(2,2,2)
        plt.hist([ item[val_names['max']] for item in data ])
        plt.title('Histogram of column "max"')
        plt.xlabel('values')
        plt.ylabel('frequency')

        plt.subplot(2,2,3)
        plt.hist([ item[val_names['min']] for item in data ])
        plt.title('Histogram of column "min"')
        plt.xlabel('values')
        plt.ylabel('frequency')

        plt.show()
    else:
        values = [ item[val_names[value]] for item in data ]
        plt.hist(values)
        plt.title('Histogram of column \"{}\" in file {}'.format(value, filename))
        plt.xlabel('values')
        plt.ylabel('frequency')
        plt.show()

def clean_and_sort_mutants(data, threshold_bio=0.2):
    
    data = [ row for row in data if row[1] == row[3] == row[5] == 'optimal' and row[2] > threshold_bio ]
    
    data.sort(key=lambda x: x[6]+(x[4]-x[6])/2)

    return data


if __name__ == '__main__':

    filenames = [
        'results_question_1_acetate.dat',
        'results_question_1_d-lactate.dat',
        'results_question_1_ethanol.dat',
        'results_question_1_l-lactate.dat',
        'results_question_1_mutial_lactate.dat',
        'results_question_1_succinate.dat']
    
    # Exchange reaction     Metabolite name
    # =====================================
    # EX_ac(e)              acetate
    # EX_lac-D(e)           D-Lactate
    # EX_lac-L(e)           L-Lactate
    # EX_succ(e)            Succinate
    # EX_etoh(e)            Ethanol
    exchanges = {
            'acetate': 'EX_ac(e)',
            'd-lactate': 'EX_lac-D(e)',
            'l-lactate': 'EX_lac-L(e)',
            'succinate': 'EX_succ(e)',
            'ethanol': 'EX_etoh(e)'}

    print("Load modified model...")
    model = cobra.io.load_json_model('result_q1_modified_model.json')

    # Create mutual objective for both lactate
    mutual_lactate = model.problem.Objective(
            model.reactions.get_by_id("EX_lac-D(e)").flux_expression +
            model.reactions.get_by_id("EX_lac-L(e)").flux_expression)
    exchanges['mutial_lactate'] = mutual_lactate

    #inspect_file(filenames[0])

    name, data = load_pickle_file(filenames[0])
    candidates = clean_and_sort_mutants(data)

    num_export_candidates = 2
    for key, exc in exchanges.items():
        csv_export_data = [['candidate', 'genes', 'max growth', 'max flux', 'min flux']]
        for idx in range(1,num_export_candidates+1):
            fig_title = "Mutant optimized for {}, candidate #{}. \nDeactivated genes: ".format(key, idx)
            gene_ids = candidates[-idx][0]
            genes_text = gene_ids[0]
            for idx in range(1,len(gene_ids)):
                genes_text += ', ' + gene_id
            fig_title += genes_text
            pdf_title = "result_q1_{}_{}".format(key, idx)
            draw_maxmin_by_growth(model, candidates[-idx][0], exc, 100, title=fig_title, to_pdf=True, pdf_title=pdf_title)
            csv_export_data.append([idx, genes_text, candidates[-idx][2], candidates[-idx][4], candidates[-idx][6]])
        
        # export csv for current exchange reaction
        with open('result_q1_{}.csv'.format(key), 'w') as csv_f:
            wr = csv.writer(csv_f, quoting=csv.QUOTE_ALL)
            wr.writerows(csv_export_data)
        
    # Export data of wild type
    csv_export_data = [['ex. reac.', 'max flux', 'min flux']]
    for key, exc in exchanges.items():
        fig_title = "Results of 'wild type' bacteria, \noptimized for {}".format(key)
        pdf_title = "result_q1_{}_wild_type".format(key)
        draw_maxmin_by_growth(model, [], exc, 100, title=fig_title, to_pdf=True, pdf_title=pdf_title)
        #csv_export_data.append([idx, genes_text, candidates[-idx][2], candidates[-idx][4], candidates[-idx][6]])
    
    # export csv for current exchange reaction
    with open('result_q1_wild_type.csv', 'w') as csv_f:
        wr = csv.writer(csv_f, quoting=csv.QUOTE_ALL)
        wr.writerows(csv_export_data)

    

    #print(candidates)

    #inspect_data(candidates)
    #draw_histogram_of_data(candidates)

    
    #draw_histogram_of_file(filenames[0], value='bio')
    #draw_histogram_of_file(filenames[0], value='max')
    #draw_histogram_of_file(filenames[0], value='min')
    
    print("END")
