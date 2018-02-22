
import cobra.test
import os
from os.path import join
import itertools




model = cobra.io.load_matlab_model(join("Model_iJO1366.mat"))

# Create list of knock out gene pairs
#gene_ids = [ gene.id for gene in model.genes ]
#gene_pairs = []
#for m in range(0, len(gene_ids)):
#    for n in range(m+1, len(gene_ids)):
#        gene_pairs.append([gene_ids[m], gene_ids[n])
gene_pairs = itertools.combinations(model.genes, 2)


def calc_maxmin_of_metabolite(gene_pairs, metabolite_id):
    # list of results, each element consists of
    # (0) maximum biomass growth
    # (1) minimum flux of metabolite when biomass is fixed to maximum
    # (2) maximum flux of metabolite when biomass is fixed to maximum
    results = []
    for gene_pair in gene_pairs:
        with model:
            #Knock out genes in current pair
            for gene in gene_pair:
                gene.knock_out()
    
            # Opitmize for biomass
            sol = model.optimize()
            opt_bio = sol.objective_value

            if opt_bio > 0:
                # Set minimum value for biomass production
                reaction = model.reactions.get_by_id('Ec_biomass_iJO1366_core_53p95M')
                reaction.lower_bound = opt_bio

                model.objective = metabolite_id
                # Optimize for maximum flux of metabolite
                sol = model.optimize()
                max_flux = sol.objective_value
                # Optimize for minimum flux of metabolite
                model.optimize.direction = 'min'
                sol = model.optimize()
                min_flux = sol.opjective_value

                results.append([gene_pair, opt_bio, min_flux, max_flux])

    return results





