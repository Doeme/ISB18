
import cobra.test
import os
from os.path import join
import itertools
import pickle
from pickle import UnpicklingError
import sys
import math
import optlang.interface
from multiprocessing.dummy import Pool as ThreadPool
from multiprocessing import Pool as ProcessPool
from multiprocessing import Value
from threading import Lock
from functools import partial
import time
from q1_simplify_model import simplify_model_with_fva, print_model_properties

# Notes:
# 
# How further boost the script?
#    > in cobrapy implementation of double_gene_deletion() a sub(sub...)call of
#      _gene_deletion() calls find_gene_knockout_reactions()
#      (defined in https://github.com/opencobra/cobrapy/blob/devel/cobra/flux_analysis/deletion.py)
#      to find all reactions which are turned off by knocking out a single gene. 
#      This function is defined in 
#      https://github.com/opencobra/cobrapy/blob/devel/cobra/manipulation/delete.py
#      and shuld boost the script additionally, I don't know why, but it would be
#      worth a try

def _init_worker(model, bio_threshold, exchange_threshold, exchange_id, pc):
    global _model
    global _bio_threshold
    global _exchange_threshold
    global _exchange_id
    global _pc

    _model = model
    _bio_threshold = bio_threshold
    _exchange_threshold = exchange_threshold
    _exchange_id = exchange_id
    _pc = pc

def _worker_thread(gene_pair):

    global _model
    global _bio_threshold
    global _exchange_threshold
    global _exchange_id
    global _pc

    # Worker thread returns a list with the folling content:
    # (0) list with knocked out genes
    # (1) solver status after max biomass optimization
    # (2) maximum biomass growth
    # (3) solver status after max flux optimization
    # (4) maximum flux of metabolite when biomass is fixed to maximum
    # (5) solver status after min flux optimization
    # (6) minimum flux of metabolite when biomass is fixed to maximum

    with _model:

        # Configure solver timeout (milliseconds)
        _model.solver.configuration.timeout = 30 * 1000
    
        # Knock out genes in current pair
        for gene_id in gene_pair:
            if not gene_id == '':
                cur_gene = _model.genes.get_by_id(gene_id)
                cur_gene.knock_out()
    
        # Opitmize for biomass
        opt_bio = _model.slim_optimize()
        opt_bio_status = _model.solver.status
    
        # Proceed if optimization solver was successful
        if opt_bio_status == optlang.interface.OPTIMAL:
            if opt_bio > _bio_threshold:
    
                # Set minimum value for biomass production
                #reaction = _model.reactions.get_by_id('Ec_biomass_iJO1366_core_53p95M')
                reaction = model.reactions.get_by_id('Biomass_Ecoli_core')
                reaction.lower_bound = opt_bio
    
                _model.objective = _exchange_id
                # Optimize for maximum flux of exchange reaction
                max_flux = _model.slim_optimize()
                max_flux_status = _model.solver.status
                if max_flux_status == optlang.interface.OPTIMAL:
                    if max_flux > _exchange_threshold:
                        # Optimize for minimum flux of exchange reaction
                        _model.objective_direction = 'min'
                        min_flux = model.slim_optimize()
                        min_flux_status = _model.solver.status
                    else:
                        max_flux_status = 'threshold'
                        min_flux = 'None'
                        min_flux_status = 'None'
                else:
                    min_flux = 'None'
                    min_flux_status = 'None'
            else:
                opt_bio_status = 'threshold'
                max_flux = 'None'
                max_flux_status = 'None'
                min_flux = 'None'
                min_flux_status = 'None'
        else:
            max_flux = 'None'
            max_flux_status = 'None'
            min_flux = 'None'
            min_flux_status = 'None'

    with _pc.get_lock():
        _pc.value += 1

    return [gene_pair,
            opt_bio_status,
            opt_bio,
            max_flux_status,
            max_flux,
            min_flux_status,
            min_flux]


def calc_maxmin_of_exchange_reaction(model, gene_pairs, exchange_id, bio_threshold, exchange_threshold, simulation_id):

    print("Generate max/min values of exchange reaction")

    # Create progress counter
    # 'I' determines the datatype, see
    # https://docs.python.org/2/library/multiprocessing.html#multiprocessing.Value
    # https://docs.python.org/2/library/array.html#module-array
    pc = Value('I', 0)

    with model:
#        pool = ThreadPool(8, initializer=_init_worker, initargs=(model, threshold, exchange_id,))
        pool = ProcessPool(processes=4, initializer=_init_worker, initargs=(model, bio_threshold, exchange_threshold, exchange_id, pc))
    
        filename = "results_question_1_{}.dat".format(simulation_id)
        with open(filename, 'wb') as results_file:
            print("   Write results to \"{}\"".format(filename))
            pickle.dump(simulation_id, results_file)
    
            results = pool.map_async(
                    _worker_thread,
                    gene_pairs)
    
            # Print progress
            gene_pairs_len = len(gene_pairs)
            while not results.ready():
                with pc.get_lock():
                    sys.stdout.write('   Genes: %s/%s (%s%s)\r' % 
                        (pc.value, gene_pairs_len, int(math.ceil(pc.value/gene_pairs_len*100)), '%'))
                time.sleep(1/30)
            sys.stdout.write('   Genes: %s/%s (%s%s)\r' % 
                (pc.value, gene_pairs_len, int(math.ceil(pc.value/gene_pairs_len*100)), '%'))
            sys.stdout.write('\n')
    
            # Check if all threads returned successfully
            if results.successful():
                print("   Threads returned without errors")
            else:
                print("   Warning: Threads returned with errors!")
    
            # Save results to file
            pickle.dump(results.get(), results_file)
    
    return results.get()


if __name__ == '__main__':

    # Configuration
    # =======================================
    path_to_models = "./../"
    model_filename = "Model_iJO1366.mat"
    gene_knock_outs = 1
    include_wild_type_model = True
    simplify_model = False
    # Only mutants with a greater bio growth will be considered
    # (in percent, relative to max of wild type)
    bio_threshold = 20
    # Only mutants with a greater exchange flow will be considered
    # (in percent, relative to max of wild type)
    exchange_threshold = 100
    # Used insead of the thresholds above if one equals zero
    # (set to zero, if zero should be used)
    approx_zero_threshold = 1e-9
    # =======================================
    
    print("Load matlab model...")
    #model = cobra.io.load_matlab_model(join(path_to_models, model_filename))
    model = cobra.test.create_test_model("textbook")
    
    # Set input fluxes of o2 and gcl to maximum
#    EX_o2 = model.reactions.get_by_id('EX_o2(e)')
#    EX_o2.lower_bound = -1000
#    EX_glc = model.reactions.get_by_id('EX_glc(e)')
#    EX_glc.lower_bound = -1000
    
    if simplify_model:
        print("Model before reaction deletion:")
        print("===============================")
        print_model_properties(model)
        print("\n")
        print("Simplifying model...")
        del_reactions, fva_results = simplify_model_with_fva(model, threshold=1e-11)
        print("Model after reaction deletion:")
        print("==============================")
        print_model_properties(model)
    
    # Create list of genes
    gene_ids = [gene.id for gene in model.genes]
    if include_wild_type_model:
        gene_ids.append('')

#    print(gene_ids)
    # Create list of knock out gene pairs
    gene_pairs = list(itertools.combinations(gene_ids, gene_knock_outs))
#    print(gene_pairs)

#    if include_lower_order_knock_outs:
#        for n in range(1,gene_knock_outs):
#            gene_pairs.append(list(itertools.combinations( [ gene.id for gene in model.genes ], n )))
    
    # Shorten the list a bit during debugging
    #gene_pairs = gene_pairs[:100]
    
    # Exchange reaction     Metabolite name
    # =====================================
    # EX_ac(e)              acetate
    # EX_lac-D(e)           D-Lactate
    # EX_lac-L(e)           L-Lactate
    # EX_succ(e)            Succinate
    # EX_etoh(e)            Ethanol
    #exchanges = {
    #        'acetate': 'EX_ac(e)', 
    #        'd-lactate': 'EX_lac-D(e)', 
    #        'l-lactate': 'EX_lac-L(e)', 
    #        'succinate': 'EX_succ(e)', 
    #        'ethanol': 'EX_etoh(e)'}
    exchanges = {
            'acetate': 'EX_ac_e', 
            'd-lactate': 'EX_lac__D_e', 
            'succinate': 'EX_succ_e', 
            'ethanol': 'EX_etoh_e'}
    
    
    # Validate the exchange ids
    for key, exc in exchanges.items():
        try:
            model.reactions.get_by_id(exc)
        except KeyError:
            sys.exit("Error: Exchange id \"{}\" is not valid".format(exc))
    
    # Create mutual objective for both lactate
#    mutual_lactate = model.problem.Objective(
#            model.reactions.get_by_id("EX_lac-D(e)").flux_expression + 
#            model.reactions.get_by_id("EX_lac-L(e)").flux_expression)
#    exchanges['mutial_lactate'] = mutual_lactate

    # Save modified model
    cobra.io.save_json_model(model, 'result_q1_modified_model.json')

    # Calculate maximum bio growth of wild type model
    with model:
        bio_ref = model.slim_optimize()
        if not model.solver.status == optlang.interface.OPTIMAL:
            sys.exit("Error: Solver was not successful when calculating the maximum bio growth of the wild type model")

    # Calculate maximum fluxes of all exchange reactions
    exchanges_ref = {}
    for key, exc in exchanges.items():
        with model:
            # Set minimum value for biomass production
            #reaction = model.reactions.get_by_id('Ec_biomass_iJO1366_core_53p95M')
            reaction = model.reactions.get_by_id('Biomass_Ecoli_core')
            reaction.lower_bound = bio_ref

            model.objective = exc
            exchanges_ref[key] = model.slim_optimize()
            if not model.solver.status == optlang.interface.OPTIMAL:
                sys.exit("Error: Solver was not successful when calculating the maximum flux of \"{}\" of the wild type model".format(key))
    
    print("Calculated reference fluxes:")
    print("  bio: {}".format(bio_ref))
    for key, ref in exchanges_ref.items():
        print("  {}: {}".format(key, ref))

    # Calculate bio growth threshold value used to kick out bad candidates
    if bio_threshold > 0:
        bio_threshold_value = bio_ref*bio_threshold/100
    else:
        if approx_zero_threshold > 0:
            bio_threshold_value = approx_zero_threshold
        else:
            bio_threshold_value = 0
    print("Used threshold values:")
    print("  bio: {}".format(bio_threshold_value))

    # Calculate exchange flux threshold value used to kick out bad candidates
    exchange_threshold_value = {}
    for key, ref in exchanges_ref.items():
        if exchange_threshold > 0:
            exchange_threshold_value[key] = exchanges_ref[key]*exchange_threshold/100
        else:
            if approx_zero_threshold > 0:
                exchange_threshold_value[key] = approx_zero_threshold
            else:
                exchange_threshold_value[key] = 0
        print("  {}: {}".format(key, exchange_threshold_value[key]))

    for key, exc in exchanges.items():
        print("Simulate mutants for exchange reaction \"{}\"".format(key))
        calc_maxmin_of_exchange_reaction(model, gene_pairs, exc, bio_threshold_value, exchange_threshold_value[key], key)
    
    print("END")
    
    
    
