# Note: It would be a good idea to create a more user-friendly version 
# of add_reaction (a la Matlab Cobra) - if it doesn't exist already 
import cobra.test 
from cobra import Model, Reaction, Metabolite 
import numpy;

import matplotlib.pyplot as plt

model = Model('hd5_3') 
 
A = Metabolite('A') 
B = Metabolite('B') 
C = Metabolite('C') 
 
# Adding "uptake" reaction b1 : -> A 
reaction = Reaction('b1') 
reaction.add_metabolites({A: -1.0}) # This enables influx of A into the network 
reaction.lower_bound = -1000. 
model.add_reaction(reaction) 
 
# Add V1 : A -> B to the network 
reaction = Reaction('V1') 
reaction.add_metabolites({A: -1.0, B: 1.0}) 
model.add_reaction(reaction) # default lower bounds=0 and upper bounds=1000 (forward reaction) 
 
# Add V2 : A -> C 
reaction = Reaction('V2') 
reaction.add_metabolites({A: -1.0, C: 1.0}) 
model.add_reaction(reaction) 
 
# Add V3 : C -> A 
reaction = Reaction('V3') 
reaction.add_metabolites({C: -1.0, A: 1.0}) 
model.add_reaction(reaction) 
 
# Add V4 : C -> B 
reaction = Reaction('V4') 
reaction.add_metabolites({C: -1.0, B: 1.0}) 
model.add_reaction(reaction) 
 
# Adding "secretion" reaction b2 : B -> 
reaction = Reaction('b2') 
reaction.add_metabolites({B: -1.0}) 
model.add_reaction(reaction) 
 
# Adding "secretion" reaction b3 : C -> 
reaction = Reaction('b3') 
reaction.add_metabolites({C: -1.0}) 
model.add_reaction(reaction) 


def print_model():
	print('%i reactions' % len(model.reactions)) 
	print('%i metabolites' % len(model.metabolites)) 
	print('%i genes' % len(model.genes)) 
	
	print("Reactions") 
	print("---------") 
	for x in model.reactions: 
		print("%s : %s" % (x.id, x.reaction)) 

print_model();
	
def optimize():
	sol = model.optimize();
	print("Fluxes: {}".format(sol.fluxes));
	print("Objective Value: {}".format(sol.objective_value));
	print("Solver Status: {}".format(sol.status));

# a)
print()
print("a)");
optimize();

# b)
model.objective='b2'
sol = model.optimize();
print();
print("b)");
print("Optimizing for flux b2");
optimize();

# c)
print();
print("c) ");
model.reactions.get_by_id("b1").lower_bound=-14
optimize();
print();
model.summary();

# d)
print();
print("d) ");
model.objective='b2'
fva_result = cobra.flux_analysis.flux_variability_analysis(model);
import pandas
print(pandas.DataFrame.from_dict(fva_result).T.round(5));

results=[];
points=numpy.linspace(0,14.0,20)

b3=model.reactions.get_by_id("b3"); 

for lb in points:
	b3.lower_bound=lb;
	b3.upper_bound=lb;
	sol = model.optimize();
	results.append(sol.objective_value)

#import matplotlib.pyplot as plt
#for t in pairs:
#	print("{}, {}".format(t[0],t[1]));

#Plot list
plt.plot(points, results, 'bo', points, results, 'k')
plt.title('Varying flux in B3')
plt.xlabel('flux in B3')
plt.ylabel('opt. flux in B2')
plt.show()

# What effect does the "mutation" have on the flow and target values? Slightly twisted
# Additional Question: Was there any reason why the Population Algorithm finds out
# This solution before (ie before blocking V1)?
# f)
print();
model.reactions.get_by_id("V1").lower_bound = 0. 
model.reactions.get_by_id("V1").upper_bound = 0. 
b3.upper_bound=1000;
b3.lower_bound=0;
optimize();

# g)

print("Adding new reaction");
D = Metabolite('D') 

# Add V5 : 2B -> D 
V5 = Reaction('V5') 
V5.add_metabolites({B: -2.0, D: 1.0}) 
model.add_reaction(V5)

# Add B4
B4 = Reaction('B4') 
B4.add_metabolites({D: -1.0}) 
model.add_reaction(B4)

print();
print_model();

results=[];
points=numpy.linspace(0,7.0,20)

for lb in points:
	V5.lower_bound=lb;
	V5.upper_bound=lb;
	sol = model.optimize();
	results.append(sol.objective_value)

#for pair in pairs:
#	print("{}, {}".format(pair[0],pair[1]));

#Plot list
plt.plot(points, results, 'bo', points, results, 'k')
plt.title('Varying flux in B5')
plt.xlabel('flux in B5')
plt.ylabel('opt. flux in B2')
plt.show()
