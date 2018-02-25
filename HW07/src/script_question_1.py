
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

class test():
    def __call__(self, item):
        print("called")

def _init_worker(model, threshold, exchange_id, pc):
    global _model
    global _threshold
    global _exchange_id
    global _pc

    _model = model
    _threshold = threshold
    _exchange_id = exchange_id
    _pc = pc

def _worker_thread(gene_pair):

    global _model
    global _threshold
    global _exchange_id
    global _pc

    # Worker thread returns a list with the folling content:
    # (0) list with knocked out genes
    # (1) solver status after max biomass optimization
    # (2) maximum biomass growth
    # (3) solver status after max flux optimization
    # (4) minimum flux of metabolite when biomass is fixed to maximum
    # (5) solver status after min flux optimization
    # (6) maximum flux of metabolite when biomass is fixed to maximum

    # Configure solver timeout (milliseconds)
    _model.solver.configuration.timeout = 30 * 1000

    # Knock out genes in current pair
    for gene in gene_pair:
        gene.knock_out()

    # Opitmize for biomass
    opt_bio = _model.slim_optimize()
    opt_bio_status = _model.solver.status

    # Proceed if optimization solver was successful
    if opt_bio_status == optlang.interface.OPTIMAL:
        if opt_bio > _threshold:

            # Set minimum value for biomass production
            reaction = _model.reactions.get_by_id('Ec_biomass_iJO1366_core_53p95M')
            reaction.lower_bound = opt_bio

            _model.objective = _exchange_id
            # Optimize for maximum flux of exchange reaction
            max_flux = _model.slim_optimize()
            max_flux_status = _model.solver.status
            # Optimize for minimum flux of exchange reaction
            _model.objective_direction = 'min'
            min_flux = model.slim_optimize()
            min_flux_status = _model.solver.status
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


def calc_maxmin_of_exchange_reaction(model, gene_pairs, exchange_id, threshold, simulation_id):

    print("Generate max/min values of exchange reaction")

    # Create progress counter
    # 'I' determines the datatype, see
    # https://docs.python.org/2/library/multiprocessing.html#multiprocessing.Value
    # https://docs.python.org/2/library/array.html#module-array
    pc = Value('I', 0)

    with model:
#        pool = ThreadPool(8, initializer=_init_worker, initargs=(model, threshold, exchange_id,))
        pool = ProcessPool(processes=4, initializer=_init_worker, initargs=(model, threshold, exchange_id, pc))
    
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



path_to_models = "./../"

print("Load matlab model...")
model = cobra.io.load_matlab_model(join(path_to_models, "Model_iJO1366.mat"))

# Create list of knock out gene pairs
gene_pairs = list(itertools.combinations(model.genes, 2))

# Shorten the list a bit during debugging
gene_pairs = gene_pairs[:100]

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


# Validate the exchange ids
for key, exc in exchanges.items():
    try:
        model.reactions.get_by_id(exc)
    except KeyError:
        sys.exit("Error: Exchange id \"{}\" is not valid".format(exc))

# Create mutual objective for both lactate
mutual_lactate = model.problem.Objective(
        model.reactions.get_by_id("EX_lac-D(e)").flux_expression + 
        model.reactions.get_by_id("EX_lac-L(e)").flux_expression)
exchanges['mutial_lactate'] = mutual_lactate


for key, exc in exchanges.items():
    print("Load exchange reaction \"{}\"".format(key))
    calc_maxmin_of_exchange_reaction(model, gene_pairs, exc, 0, key)

print("END")



