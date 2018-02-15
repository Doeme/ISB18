#!/usr/bin/python3
import cobra;
import cobra.test;

import sys;
epsilon=sys.float_info.epsilon


def generate_escher(filename, sol):
	f=open(filename,"w+");
	f.write(sol.fluxes.to_csv());
	f.close();

#ecoli=cobra.io.load_matlab_model("../EColi.mat")
ecoli=cobra.test.create_test_model("ecoli");

print()
print("2.a)")
ex_O2=ecoli.reactions.get_by_id("EX_o2_e");
ex_O2.lower_bound=-1000;
sol=ecoli.optimize()
generate_escher("2_a.csv",sol);
o2=sol.fluxes['EX_o2_e']
glc=sol.fluxes['EX_glc_e'];
print("i) 02/glc = {}".format(o2/glc))
print("Biomass Production (Objective value): {}".format(sol.objective_value))
sol_aerobic=sol;

print()
print("2.b)")
ex_O2.lower_bound=0;
sol=ecoli.optimize()
generate_escher("2_b.csv",sol);
o2=sol.fluxes['EX_o2_e']
glc=sol.fluxes['EX_glc_e'];
print("i) 02/glc = {}".format(o2/glc))
print("ii) Biomass Production (Objective value): {}".format(sol.objective_value))
exchanges=[r.id for r in ecoli.exchanges];
exchange_fluxes=sol.fluxes[exchanges];
secretions=exchange_fluxes.select(lambda x: exchange_fluxes[x]>0);
print("iii) Secretions: ");
print(secretions.to_string());

sol_anaerobic=sol;

import cobra.flux_analysis;

print();
print("3.a) Forward reactions:");
fva = cobra.flux_analysis.flux_variability_analysis;
res = fva(ecoli, ecoli.reactions);
fwd_reactions = res.select(lambda x: res['minimum'][x] >= epsilon);
print(fwd_reactions.to_string());

print();
print("3.b) Forced secretions:");
exchanges=[r.id for r in ecoli.exchanges];
exchange_fluxes = res.select(lambda x: x in exchanges);
forced_secretion = exchange_fluxes.select(lambda x: res['minimum'][x] > epsilon);
print(forced_secretion.to_string());

print();
print("3.c) Blocked reactions:");
blocked = res.select(lambda x: (abs(res['minimum'][x]) < epsilon and abs(res['maximum'][x]) < epsilon) );
print(blocked.to_string());

print();
print("3.d) Reactions with large variance. This can come from cycles in the graph");
large_variance = [a for a,b,c, in res.itertuples() if abs(b-c)>10];
large_variance = res.select(lambda x: x in large_variance);
print(large_variance.to_string());

#from matplotlib import *;
import matplotlib.pyplot as plt

print();
print("3.e) ");
f=open("barplot.tex","w+");
for a,b,c in res.itertuples():
	b=min(35,b);
	c=max(-25,c);
	print("\\verb|{0}| & \\barplot{{{1:.4f}}}{{{2:.4f}}}\\\\".format(a,b,c),file=f)
f.close();
