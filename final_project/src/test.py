import dmmm;
import cobra;

l=cobra.io.read_sbml_model("../models/lactobacillus.sbml")
l.name="Lactobacillus Brevis";
ll=dmmm.Organism(l,1,0.01)
ll.exchange_filter=lambda x: x.id!="added_biomass_sink";

y=cobra.io.read_sbml_model("../models/yeast.sbml")
y.name="Sacharomyces Cerevisiæ";
yy=dmmm.Organism(y,1000,0.01)

co=[ll,yy]

d=dmmm.DMMM(co)


vol=50;
mol_vol=55346;
perc_glc=0.10;

d.medium.metabolites["EX_h2o_e_"]=(1-perc_glc)*vol*mol_vol*1000;
d.medium.metabolites["EX_glc_e_"]=(perc_glc)*vol*mol_vol*1000;
d.medium.metabolites["EX_nh4_e_"]=(perc_glc)*vol*mol_vol*1000*0.1;
d.medium.metabolites["EX_so4_e_"]=vol*mol_vol*1000*0.0001;
d.medium.metabolites["EX_pi_e_"]=vol*mol_vol*1000*0.01;
d.medium.metabolites["EX_Nbfortyr_e_"]=vol*mol_vol*1000*0.01;
d.medium.metabolites["EX_o2_e_"]=vol*mol_vol*1000*0.01;

trace_elements=['EX_arab_L_e_', 'EX_cys_L_e_', 'EX_glu_L_e_', 'EX_ile_L_e_', 'EX_leu_L_e_', 'EX_phe_L_e_', 'EX_val_L_e_']
for id in trace_elements:
	d.medium.metabolites[id]=vol*mol_vol*1000*0.001;

r=d.simulate(3*31*24);
import matplotlib.pyplot as pp;

def plot_met(r, met):
	pp.plot(r.get_times(),r.get_metabolite(met));

def plot_pop(r, idx):
	pp.plot(r.get_times(),r.get_population(idx));

import IPython; IPython.embed()
