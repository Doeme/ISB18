import cobra
from cobra.flux_analysis import flux_variability_analysis
from os.path import join
import matplotlib.pyplot as plt
import math


def pandas_to_list(pandas):
    return [ list(t) for t in pandas.itertuples() ]

def pandas_to_dict(pandas):
    return { t[0]: { 'max': t[1], 'min': t[2] } for t in pandas.itertuples() }

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

def plot_fva_with_threshold(results, threshold):
    
    res_lst = pandas_to_list(results)
    maxima = np.abs([ res[1] for res in res_lst ])
    minima = np.abs([ res[2] for res in res_lst ])

    # Remove all zero elements
    maxima = np.array([ maxima[idx] for idx in maxima.nonzero() ])
    maxima_log = np.log10(maxima)
    minima = np.array([ minima[idx] for idx in minima.nonzero() ])
    minima_log = np.log10(minima)

    # Prepare axis
#    max_end = np.ceil(np.log10(maxima.max()))
#    max_start = np.floor(np.log10(maxima.min()))
#    min_end = np.ceil(np.log10(minima.max()))
#    min_start = np.floor(np.log10(minima.min()))

#    max_axis = np.logspace(max_start, max_end, 1000)
#    min_axis = np.logspace(min_start, min_end, 1000)

    plt.subplot(2,1,1)
    plt.hist(maxima_log, bins=1000)
    plt.title('maxima')

    plt.subplot(2,1,2)
    plt.hist(minima_log, bins=1000)
    plt.title('minima')

    plt.show()



def find_unused_reactions(model, threshold=1e-9):
    with model:
        presults = flux_variability_analysis(model, model.reactions, fraction_of_optimum=0)
        results = pandas_to_list(presults)

        reaction_ids = []
        for res in results:
            if math.fabs(res[1]) < threshold and math.fabs(res[2]) < threshold:
                reaction_ids.append(res[0])
                #print("max: {}, min: {}".format(res[1], res[2]))

    return reaction_ids, presults

def simplify_model_with_fva(model, threshold=1e-9):

    reaction_ids, results = find_unused_reactions(model, threshold)
    reactions = [ model.reactions.get_by_id(r_id) for r_id in reaction_ids ]
    model.remove_reactions(reactions, remove_orphans=True)

    return reactions, results

def plot_fva_results_noise(model, reaction_ids, repetitions=10):
    
    results = []
    for i in range(0, repetitions):
        with model:
            reactions = [ model.reactions.get_by_id(ex_id) for ex_id in reaction_ids ]
            results.append(pandas_to_dict(flux_variability_analysis(model, reactions)))

    values = { r_id: {'max': [ res[r_id]['max'] for res in results ], 'min': [ res[r_id]['min'] for res in results ] } for r_id in reaction_ids }

    plt.subplot(2,2,1)
    plt.hist(values['EX_succ(e)']['max'])
    plt.title('EX_succ(e) maxima')

    plt.subplot(2,2,3)
    plt.hist(values['EX_succ(e)']['min'])
    plt.title('EX_succ(e) minima')

    plt.subplot(2,2,2)
    plt.hist(values['EX_etoh(e)']['max'])
    plt.title('EX_etoh(e) maxima')

    plt.subplot(2,2,4)
    plt.hist(values['EX_etoh(e)']['min'])
    plt.title('EX_etoh(e) minima')
    
    plt.show()

if __name__ == "__main__":


    path_to_models = './../'

    print("Load model...")
    model = cobra.io.load_json_model(join(path_to_models, 'ecolicore.json'))

    # Set input fluxes of o2 and gcl to maximum
    EX_o2 = model.reactions.get_by_id('EX_o2(e)')
    EX_o2.lower_bound = -1000
    EX_glc = model.reactions.get_by_id('EX_glc(e)')
    EX_glc.lower_bound = -1000

    exchanges = {
        'acetate': 'EX_ac(e)',
        'd-lactate': 'EX_lac-D(e)',
        'l-lactate': 'EX_lac-L(e)',
        'succinate': 'EX_succ(e)',
        'ethanol': 'EX_etoh(e)'}

#    plot_fva_results_noise(model, list(exchanges.values()) ) 

    reaction_ids, results = find_unused_reactions(model, threshold=1e-11)

    d_results = pandas_to_dict(results)

    print("Results of exchanges:")
    print("Exchange\tmax\tmin")
    print("==================")
    for key, ex_id in exchanges.items():
        print("{}\t{}\t{}".format(ex_id, d_results[ex_id]['max'], d_results[ex_id]['min']))

    del_exchanges = []
    for key, ex in exchanges.items():
        try:
            reaction_ids.index(ex)
            del_exchanges.append(ex)
        except ValueError:
            pass

#    print("Exchange ID\tStatus")
#    print("===================")
#    for idx in range(0, len(exchanges)):
#        print("{}\t{}".format(exchanges

    print("Deleted reactions:")
    print(del_exchanges)







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





