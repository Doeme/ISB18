
import cobra
import os
from os.path import join
import optlang.interface
import numpy as np
import matplotlib.pyplot as plt





def draw_maxmin_by_growth(model, gene_pair, exchange_id, points=10):


    with model:

        # Configure solver timeout (milliseconds)
        model.solver.configuration.timeout = 30 * 1000

        # Knock out genes
        for gene in gene_pair:
            gene.knock_out()

        # Optimize for biomass
        bio_max = model.slim_optimize()
        
        if model.solver.status == optlang.interface.OPTIMAL:
            
            bio_points = points
            bio_delta = bio_max/(bio_points-1)
            bio_axis = np.arange(0, bio_max+bio_delta/2, bio_delta)
            bio_axis[bio_points-1] = bio_max

            print("bio_axis: {}".format(bio_axis))

            min_fluxes = []
            max_fluxes = []

            for bio_val in bio_axis:
                # Set min/max value for biomass production
                reaction = model.reactions.get_by_id('Ec_biomass_iJO1366_core_53p95M')
                reaction.lower_bound = bio_val
                reaction.upper_bound = bio_val

                model.objective_direction = 'max'

                # Set exchange_id as new optimization objective
                model.objective = exchange_id
                # Optimize for max flux of given exchange reaction
                max_flux = model.slim_optimize()
                if model.solver.status == optlang.interface.OPTIMAL:
                    max_fluxes.append(max_flux)
                else:
                    max_fluxes.append(float('nan'))
                    print("Warning: Could not optimize for max flux in point {}. Solver returned with status \"{}\".".format(bio_val, model.solver.status))

                # Optimize for min flux of given exchange reaction
                model.objective_direction = 'min'
                min_flux = model.slim_optimize()
                if model.solver.status == optlang.interface.OPTIMAL:
                    min_fluxes.append(min_flux)
                else:
                    min_fluxes.append(float('nan'))
                    print("Warning: Could not optimize for max flux in point {}. Solver returned with status \"{}\".".format(bio_val, model.solver.status))

        else:
            print("Error: Could not optimize for biomass. Solver returned with status \"{}\".".format(model.solver.status))

        print("bio\tmin\tmax")
        for idx in range(0,bio_points):
            print("{}\t{}\t{}".format(bio_axis[idx], min_fluxes[idx], max_fluxes[idx]))

        #y_pos = np.arange(bio_points)
        plt.plot(bio_axis, max_fluxes, 'r',
                bio_axis, min_fluxes, 'g')
        plt.fill_between(bio_axis, max_fluxes, min_fluxes, color='blue', alpha=0.5)
        plt.title('Super toller plot')
        plt.xlabel('growth')
        plt.ylabel('flux')
        plt.show()

               

path_to_models = "./../"

print("Load matlab model...")
model = cobra.io.load_matlab_model(join(path_to_models, "Model_iJO1366.mat"))

exchanges = {
        'acetate': 'EX_ac(e)',
        'd-lactate': 'EX_lac-D(e)',
        'l-lactate': 'EX_lac-L(e)',
        'succinate': 'EX_succ(e)',
        'ethanol': 'EX_etoh(e)'}
#for key, exchange_id in exchanges.items():
#    print(key)

print("Calculate values...")
draw_maxmin_by_growth(model, [], exchanges['acetate'], points=100)



