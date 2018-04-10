import dmmm;
import cobra;

l=cobra.io.read_sbml_model("../models/lactobacillus.sbml")
ll=dmmm.Organism(l,1,0.0)
ll.exchange_filter=lambda x: x.id!="added_biomass_sink";

y=cobra.io.read_sbml_model("../models/yeast.sbml")
yy=dmmm.Organism(y,1000,0.0)

co=[ll,yy]

d=dmmm.DMMM(co)


for k,v in d.medium.metabolites.items():
	d.medium.metabolites[k]=v+100;

vol=50;
mol_vol=55346;
perc_glc=0.10;

d.medium.metabolites["EX_h2o_e_"]=(1-perc_glc)*vol*mol_vol*1000;
d.medium.metabolites["EX_glc_e_"]=(perc_glc)*vol*mol_vol*1000;
d.medium.metabolites["EX_o2_e_"]=(perc_glc)*vol*mol_vol*1000*10;


l=d.simulate(14*24);

import IPython; IPython.embed()
