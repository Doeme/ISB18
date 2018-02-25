
import cobra.test
import os
from os.path import join
import itertools
import pickle
from pickle import UnpicklingError
import sys
import math
import optlang.interface


def calc_maxmin_of_exchange_reaction(model, gene_pairs, exchange_id, threshold, simulation_id):

    print("Generate max/min values of exchange reaction")

    filename = "results_question_1_{}.dat".format(simulation_id)
    with open(filename, 'wb') as results_file:
        print("   Write results to \"{}\"".format(filename))
        pickle.dump(simulation_id, results_file)

        # list of results, each element consists of
        # (0) maximum biomass growth
        # (1) minimum flux of metabolite when biomass is fixed to maximum
        # (2) maximum flux of metabolite when biomass is fixed to maximum
        results = []
        gene_pairs_len = len(gene_pairs)
        idx = 0
        hits = 0
        sys.stdout.write('   Genes: %s/%s (%s%s), hits: %s\r' % 
                (idx, gene_pairs_len, int(math.ceil(idx/gene_pairs_len*100)), '%', hits))
        for gene_pair in gene_pairs:
            with model:
                #Knock out genes in current pair
                for gene in gene_pair:
                    gene.knock_out()
        
                # Opitmize for biomass
                opt_bio = model.slim_optimize()
                opt_bio_status = model.solver.status
    
                # Proceed if optimization solver was successful
                if opt_bio_status == optlang.interface.OPTIMAL:
                    if opt_bio > threshold:
                        hits += 1
                        # Set minimum value for biomass production
                        reaction = model.reactions.get_by_id('Ec_biomass_iJO1366_core_53p95M')
                        reaction.lower_bound = opt_bio
        
                        model.objective = exchange_id
                        # Optimize for maximum flux of exchange reaction
                        max_flux = model.slim_optimize()
                        max_flux_status = model.solver.status
                        # Optimize for minimum flux of exchange reaction
                        model.objective_direction = 'min'
                        min_flux = model.slim_optimize()
                        min_flux_status = model.solver.status
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

                # Save results to file and flush all buffers
                pickle.dump([
                    gene_pair, 
                    opt_bio_status, 
                    opt_bio, 
                    min_flux_status, 
                    min_flux, 
                    max_flux_status, 
                    max_flux], results_file)
                results_file.flush()
                os.fsync(results_file.fileno())
    
                # Update progress
                idx += 1
                sys.stdout.flush()
                sys.stdout.write('   Genes: %s/%s (%s%s), hits: %s\r' % 
                    (idx, gene_pairs_len, int(math.ceil(idx/gene_pairs_len*100)), '%', hits))

    

    sys.stdout.write('\n')
    return results



path_to_models = "./../"

print("Load matlab model...")
model = cobra.io.load_matlab_model(join(path_to_models, "Model_iJO1366.mat"))

# Create list of knock out gene pairs
#gene_ids = [ gene.id for gene in model.genes ]
#gene_pairs = []
#for m in range(0, len(gene_ids)):
#    for n in range(m+1, len(gene_ids)):
#        gene_pairs.append([gene_ids[m], gene_ids[n])
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
#exchanges = ['EX_ac(e)', 'EX_lac-D(e)', 'EX_lac-L(e)', 'EX_succ(e)', 'EX_etoh(e)']
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



