import cobra
from cobra.flux_analysis import flux_variability_analysis
from os.path import join
import matplotlib.pyplot as plt
import math


def pandas_to_list(pandas):
    return [ list(t) for t in pandas.itertuples() ]

def compare_fva_results_with():
    path_to_models = './../'
    
    print("Load model...")
    model = cobra.io.load_json_model(join(path_to_models, 'ecolicore.json'))
    
    with model:
        print("FVA with untouched model...")
        pandas_low_flux = flux_variability_analysis(model, model.reactions)
        results_low_flux = pandas_to_list(pandas_low_flux)
    
    
    with model:
        print("FVA with changed input flux...")
        # Set input fluxes of o2 and gcl to maximum
        EX_o2 = model.reactions.get_by_id('EX_o2(e)')
        EX_o2.lower_bound = -1000
        EX_glc = model.reactions.get_by_id('EX_glc(e)')
        EX_glc.lower_bound = -1000
    
        pandas_high_flux = flux_variability_analysis(model, model.reactions)
        results_high_flux = pandas_to_list(pandas_high_flux)
    
    
    diff = pandas_low_flux.values - pandas_high_flux.values
    
    plt.subplot(1, 2, 1)
    plt.hist([r[0] for r in diff])
    plt.title('max')
    
    plt.subplot(1, 2, 2)
    plt.hist([r[1] for r in diff])
    plt.title('min')
    
    plt.show()


def find_unused_reactions(model, threshold=1e-9):
    with model:
        presults = flux_variability_analysis(model, model.reactions)
        results = pandas_to_list(presults)

        reaction_ids = []
        for res in results:
            if math.fabs(res[1]) < threshold and math.fabs(res[2]) < threshold:
                reaction_ids.append(res[0])
                #print("max: {}, min: {}".format(res[1], res[2]))

    return reaction_ids, results

def simplify_model_with_fva(model, threshold=1e-9):

    reaction_ids, results = find_unused_reactions(model, threshold)
    reactions = [ model.reactions.get_by_id(r_id) for r_id in reaction_ids ]
    model.remove_reactions(reactions, remove_orphans=True)

    return reactions, results

#path_to_models = './../'
#
#print("Load model...")
#model = cobra.io.load_json_model(join(path_to_models, 'ecolicore.json'))
#
## Set input fluxes of o2 and gcl to maximum
#EX_o2 = model.reactions.get_by_id('EX_o2(e)')
#EX_o2.lower_bound = -1000
#EX_glc = model.reactions.get_by_id('EX_glc(e)')
#EX_glc.lower_bound = -1000
#
#with model:
#    results = pandas_to_list(flux_variability_analysis(model, model.reactions))
#
#plt.subplot(1,2,1)
#n, bins, patches = plt.hist([ r[1] for r in results ], histtype='stepfilled')
#plt.title('max')
##plt.xticks(bins)
#print("max histogram:")
#print("Frequency:")
#print(n)
#print("Values:")
#print(bins)
#
#plt.subplot(1,2,2)
#n, bins, patches = plt.hist([ r[2] for r in results ], bins=1000)
#plt.title('min')
##plt.xticks(bins)
#print("min histogram:")
#print("Frequency:")
#print(n)
#print("Values:")
#print(bins)
#plt.show()


def print_model_properties(model):
    print("genes:              {}".format(len(model.genes)))
    print("reactions:          {}".format(len(model.reactions)))
    print("exchange reactions: {}".format(len(model.exchanges)))
    print("metabolites:        {}".format(len(model.exchanges)))

#print("Model before reaction deletion:")
#print("===============================")
#print_model_properties(model)
#print("\n")
#print("Delete {}/{} reactions".format(len(reactions), len(model.reactions)))
#
#model.remove_reactions(reactions, remove_orphans=True)
#
#print("Model after reaction deletion:")
#print_model_properties(model)





