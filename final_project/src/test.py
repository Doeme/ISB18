import dmmm;
import cobra;

l=cobra.io.read_sbml_model("../models/lactobacillus.sbml")
ll=dmmm.Organism(l,1,0.005)
ll.exchange_filter=lambda x: x.id!="added_biomass_sink";

y=cobra.io.read_sbml_model("../models/yeast.sbml")
yy=dmmm.Organism(y,1000,0.02)

co=[ll,yy]

d=dmmm.DMMM(co)


vol=50;
mol_vol=55346;
perc_glc=0.10;

d.medium.metabolites["EX_h2o_e_"]=(1-perc_glc)*vol*mol_vol*1000;
d.medium.metabolites["EX_glc_e_"]=(perc_glc)*vol*mol_vol*1000;
d.medium.metabolites["EX_nh4_e_"]=(perc_glc)*vol*mol_vol*1000;
d.medium.metabolites["EX_so4_e_"]=vol*mol_vol*1000*0.01;
d.medium.metabolites["EX_pi_e_"]=vol*mol_vol*1000*0.01;
d.medium.metabolites["EX_Nbfortyr_e_"]=vol*mol_vol*1000*0.01;
d.medium.metabolites["EX_o2_e_"]=vol*mol_vol*1000*0.0001;


r=d.simulate(12*7*24);
import matplotlib.pyplot as pp;

def plot_met(r, met):
	pp.plot(r.get_times(),r.get_metabolite(met));

def plot_pop(r, idx):
	pp.plot(r.get_times(),r.get_population(idx));

import IPython; IPython.embed()
