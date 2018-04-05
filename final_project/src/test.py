import dmmm;
import cobra;

l=cobra.io.read_sbml_model("../Models/lactobacillus.xml")
ll=dmmm.Organism(l,1000,0)

y=cobra.io.read_sbml_model("../Models/yeast.xml")
yy=dmmm.Organism(y,100000,0)

co=[ll,yy]

d=dmmm.DMMM(co)

d.medium.metabolites["EX_glc_e_"]=1000
d.medium.metabolites["EX_o2_e_"]=50

for k,v in d.medium.metabolites.items():
	d.medium.metabolites[k]=v+1;


l=d.simulate(0.1);
print(l);
