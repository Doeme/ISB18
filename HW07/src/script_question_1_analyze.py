
import cobra
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
import optlang.interface


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
            print("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}".format(idx, data[idx][0], data[idx][1], data[idx][2], data[idx][3], data[idx][4], data[idx][5], data[idx][6]))

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
#        'results_question_1_l-lactate.dat',
#        'results_question_1_mutial_lactate.dat',
        'results_question_1_succinate.dat']
    
    # Exchange reaction     Metabolite name
    # =====================================
    # EX_ac(e)              acetate
    # EX_lac-D(e)           D-Lactate
    # EX_lac-L(e)           L-Lactate
    # EX_succ(e)            Succinate
    # EX_etoh(e)            Ethanol
#    exchanges = {
#            'acetate': 'EX_ac(e)',
#            'd-lactate': 'EX_lac-D(e)',
#            'l-lactate': 'EX_lac-L(e)',
#            'succinate': 'EX_succ(e)',
#            'ethanol': 'EX_etoh(e)'}
    exchanges = {
            'acetate': 'EX_ac_e',
            'd-lactate': 'EX_lac__D_e',
            'succinate': 'EX_succ_e',
            'ethanol': 'EX_etoh_e'}

    print("Load modified model...")
    model = cobra.io.load_json_model('result_q1_modified_model.json')
    #model = cobra.test.create_test_model("textbook")

    # Create mutual objective for both lactate
#    mutual_lactate = model.problem.Objective(
#            model.reactions.get_by_id("EX_lac-D(e)").flux_expression +
#            model.reactions.get_by_id("EX_lac-L(e)").flux_expression)
#    exchanges['mutial_lactate'] = mutual_lactate

    #inspect_file(filenames[0])

    #name, data = load_pickle_file(filenames[0])
    #inspect_data(data)
    #candidates = clean_and_sort_mutants(data)

    num_export_candidates = 2
    for filename in filenames:
        print("Open file {}".format(filename))
        key, data = load_pickle_file(filename)
        print("Sort mutants for {}".format(key))
        candidates = clean_and_sort_mutants(data)

        print("Draw figures for {}".format(key))
        csv_export_data = [['candidate', 'genes', 'max growth', 'max flux', 'min flux']]
        for idx in range(1,num_export_candidates+1):
            fig_title = "Mutant optimized for {}, candidate #{}. \nDeactivated genes: ".format(key, idx)
            gene_ids = candidates[-idx][0]
            genes_text = gene_ids[0]
            for i in range(1,len(gene_ids)):
                genes_text += ' ' + gene_ids[i]
            fig_title += genes_text
            pdf_title = "result_q1_{}_{}".format(key, idx)
            print("Save {}.pdf".format(pdf_title))
            draw_maxmin_by_growth(model, candidates[-idx][0], exchanges[key], 100, title=fig_title, to_pdf=True, pdf_title=pdf_title)
            csv_export_data.append([idx, genes_text, candidates[-idx][2], candidates[-idx][4], candidates[-idx][6]])
        
        # export csv for current exchange reaction
        csv_filename = 'result_q1_{}.csv'.format(key)
        print("Export data of {} to file {}".format(key, csv_filename))
        with open(csv_filename, 'w') as csv_f:
            wr = csv.writer(csv_f)
            wr.writerows(csv_export_data)
    
    # Calculate data of wild type
    print("Calculate export data for wild type...")
    wild_type_max_fluxes = {}
    wild_type_min_fluxes = {}
    for key, exchange_id in exchanges.items():
       with model:

        # Configure solver timeout (milliseconds)
        model.solver.configuration.timeout = 30 * 1000

        # Optimize for biomass
        wild_type_bio_max = model.slim_optimize()

        if model.solver.status == optlang.interface.OPTIMAL:

            # Set min/max value for biomass production
            #reaction = model.reactions.get_by_id('Ec_biomass_iJO1366_core_53p95M')
            reaction = model.reactions.get_by_id('Biomass_Ecoli_core')
            reaction.lower_bound = wild_type_bio_max
            reaction.upper_bound = wild_type_bio_max

            model.objective_direction = 'max'

            # Set exchange_id as new optimization objective
            model.objective = exchange_id
            # Optimize for max flux of given exchange reaction
            max_flux = model.slim_optimize()
            if model.solver.status == optlang.interface.OPTIMAL:
                wild_type_max_fluxes[key] = max_flux
            else:
                wild_type_max_fluxes[key] = float('nan')
                print("Error: Could not optimize for max of {}. Solver returned with status \"{}\".".format(key, model.solver.status))

            # Optimize for min flux of given exchange reaction
            model.objective_direction = 'min'
            min_flux = model.slim_optimize()
            if model.solver.status == optlang.interface.OPTIMAL:
                wild_type_min_fluxes[key] = min_flux
            else:
                wild_type_min_fluxes[key] = float('nan')
                print("Error: Could not optimize for max of {}. Solver returned with status \"{}\".".format(key, model.solver.status))

    # Export data of wild type
    csv_filename = 'result_q1_w'
    print("Generate figures of wild type bacteria...")
    csv_export_data = [['Metabolite', 'max flux', 'min flux']]
    for key, exc in exchanges.items():
        fig_title = "Results of 'wild type' bacteria, \noptimized for {}. \nMax bio growth={}".format(key, wild_type_bio_max)
        pdf_title = "result_q1_{}_wild_type".format(key)
        print("Save {}.pdf".format(pdf_title))
        draw_maxmin_by_growth(model, [], exc, 100, title=fig_title, to_pdf=True, pdf_title=pdf_title)
        csv_export_data.append([key, wild_type_max_fluxes[key], wild_type_min_fluxes[key]])
    
    # export csv for current exchange reaction
    print("Export data of wild type bacteria to file {}".format(csv_filename))
    with open(csv_filename, 'w') as csv_f:
        wr = csv.writer(csv_f)
        wr.writerows(csv_export_data)

    

    #print(candidates)

    #inspect_data(candidates)
    #draw_histogram_of_data(candidates)

    
    #draw_histogram_of_file(filenames[0], value='bio')
    #draw_histogram_of_file(filenames[0], value='max')
    #draw_histogram_of_file(filenames[0], value='min')
    
    print("END")
