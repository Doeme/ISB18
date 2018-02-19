
import cobra.test
import os
from os.path import join
import matplotlib.pyplot as plt
import plotly.plotly as py
import pickle
import numpy as np

# path to input files
input_file_path = "./../"
output_file_path = "./../"

#model = cobra.io.load_matlab_model(join(input_file_path, "EColi.mat"))

#sol = model.optimize()

#print(sol.fluxes.to_csv())

#with open(output_file_path + "fluxes_test.txt", 'w+') as outfile:
#   outfile.write(sol.fluxes.to_csv())

def get_ko_gene_ids_from_model(model_file, th_biomass_percent):
    
    model = cobra.io.load_matlab_model(join(model_file))

    # Get reference biomass production without gene knock outs
    sol = model.optimize()
    bio_ref = sol.objective_value

    # Threshold for biomass production
    th_biomass = bio_ref*th_biomass_percent/100

    # Get list of gene ids in model
    gene_ids = [ gene.id for gene in model.genes ]

    ko_gene_ids = []
    # Find all ko genes in gene id list
    for gene_id in gene_ids:
        with model:
            # Knock out current gene in model
            cur_gene = model.genes.get_by_id(gene_id)
            cur_gene.knock_out()
            # optimize for biomass (biomass objective is set by default)
            sol = model.optimize()
            # add to list if biomass production is lower than threshold
            if sol.objective_value < th_biomass:
                ko_gene_ids.append(gene_id)

    return (gene_ids, ko_gene_ids)

ko_gene_ids_lst = []
gene_ids_lst = []
organisms = ["T. maritima", "Synechocystsis", "E. coli", "S. elongatus", "S. cerevisiae"]

# Investigate ko genes in Thermotoga maritima
gene_ids, ko_gene_ids = get_ko_gene_ids_from_model(input_file_path + "Thermotoga.mat", 10)
ko_gene_ids_lst.append(ko_gene_ids)
gene_ids_lst.append(gene_ids)

# Investigate ko genes in Synechocystsis
gene_ids, ko_gene_ids = get_ko_gene_ids_from_model(input_file_path + "Synechocystis.mat", 10)
ko_gene_ids_lst.append(ko_gene_ids)
gene_ids_lst.append(gene_ids)

# Investigate ko genes in E. coli
gene_ids, ko_gene_ids = get_ko_gene_ids_from_model(input_file_path + "EColi.mat", 10)
ko_gene_ids_lst.append(ko_gene_ids)
gene_ids_lst.append(gene_ids)

# Investigate ko genes in S. elongatus
gene_ids, ko_gene_ids = get_ko_gene_ids_from_model(input_file_path + "Synechococcus.mat", 10)
ko_gene_ids_lst.append(ko_gene_ids)
gene_ids_lst.append(gene_ids)

# Investigate ko genes in S. cerevisiae
gene_ids, ko_gene_ids = get_ko_gene_ids_from_model(input_file_path + "Model_yeast7.mat", 10)
ko_gene_ids_lst.append(ko_gene_ids)
gene_ids_lst.append(gene_ids)


for idx in range(0,len(organisms)):
    print("Organism: " + organisms[idx])
    print("KO genes: {}".format(len(ko_gene_ids_lst[idx])))


#gene_list = []
#for lst in ko_gene_ids:
#    gene_list = gene_list + lst

#gene_list = list(set(gene_list))

#print(len(gene_list))


#with open("tmp_objects.p", 'wb') as pickle_file:
#    pickle.dump(ko_gene_ids_lst, pickle_file)
#    pickle.dump(gene_ids_lst, pickle_file)
#    pickle.dump(organisms, pickle_file)


#with open("tmp_objects.p", 'rb') as pickle_file:
#    ko_gene_ids_lst = pickle.load(pickle_file)
#    gene_ids_lst = pickle.load(pickle_file)
#    organisms = pickle.load(pickle_file)


# Doubling times of organisms in minutes
doubling_time = [2*60, 12*60, 20, 6.5*60, 105]

# Generate bar chart
ko_gene_lens = [ len(lst) for lst in ko_gene_ids_lst ]
gene_lens = [ len(lst) for lst in gene_ids_lst ]


y = []
for idx in range(0, len(ko_gene_lens)):
    y.append(ko_gene_lens[idx]/gene_lens[idx]*100)

print("Organism\tgenes\tKO genes")
print("=================================")
for idx in range(0, len(organisms)):
    print("{} \t{}\t{}\t{}\t{}".format(organisms[idx], gene_lens[idx], ko_gene_lens[idx], y[idx], doubling_time[idx]))



y_pos = np.arange(len(organisms))

plt.bar(y_pos, y, align='center', alpha=0.5)
plt.xticks(y_pos, organisms)
plt.xlabel('Organisms')
plt.ylabel('Ratio in %')
plt.title('Ratio of lethal genes per organism')

axes2 = plt.twinx()
axes2.plot(y_pos, doubling_time, 'ro')
axes2.plot(y_pos, doubling_time, 'r:')
axes2.set_ylabel('cell division time in minutes')


plt.show()




print("END")
