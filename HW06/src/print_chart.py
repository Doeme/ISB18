import pickle
import matplotlib.pyplot as plt
import numpy as np


with open("tmp_objects.p", 'rb') as pickle_file:
    ko_gene_ids_lst = pickle.load(pickle_file)
    gene_ids_lst = pickle.load(pickle_file)
    organisms = pickle.load(pickle_file)


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
