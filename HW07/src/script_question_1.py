
import cobra.test
import os
from os.path import join
import itertools
import pickle
import sys
import math


def calc_maxmin_of_exchange_reaction(model, gene_pairs, exchange_id, threshold):

    print("Generate max/min values of exchange reaction")

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
            sol = model.optimize()
            opt_bio = sol.objective_value

            if opt_bio > threshold:
                hits += 1
                # Set minimum value for biomass production
                reaction = model.reactions.get_by_id('Ec_biomass_iJO1366_core_53p95M')
                reaction.lower_bound = opt_bio

                model.objective = exchange_id
                # Optimize for maximum flux of exchange reaction
                sol = model.optimize()
                max_flux = sol.objective_value
                # Optimize for minimum flux of exchange reaction
                model.objective_direction = 'min'
                sol = model.optimize()
                min_flux = sol.objective_value

                #gene_pair = [gene.id for gene in gene_pair]

                results.append([gene_pair, opt_bio, min_flux, max_flux])

            # Update progress
            idx += 1
            sys.stdout.flush()
            sys.stdout.write('   Genes: %s/%s (%s%s), hits: %s\r' % 
                (idx, gene_pairs_len, int(math.ceil(idx/gene_pairs_len*100)), '%', hits))


    sys.stdout.write('\n')
    return results




print("Load matlab model...")
model = cobra.io.load_matlab_model(join("Model_iJO1366.mat"))

# Create list of knock out gene pairs
#gene_ids = [ gene.id for gene in model.genes ]
#gene_pairs = []
#for m in range(0, len(gene_ids)):
#    for n in range(m+1, len(gene_ids)):
#        gene_pairs.append([gene_ids[m], gene_ids[n])
gene_pairs = list(itertools.combinations(model.genes, 2))

# Shorten the list a bit during debugging
#gene_pairs = gene_pairs[:100]

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
    results = calc_maxmin_of_exchange_reaction(model, gene_pairs, exc, 0)

    filename = "results_question_1_{}.dat".format(key)
    with open(filename, 'wb') as tmp_file:
        pickle.dump(key, tmp_file)
        pickle.dump(results, tmp_file)

    print("   Saved results in file \"{}\"".format(filename))

print("END")



