import dmmm;
import cobra;

MM={'etoh':46.07, 'glc':180.16, 'co2':44.0095, 'o2':15.99943, 'lac':90.08, 'ac':59.04}; #Molar masses in g/mol

vol=50;
mol_vol=55.408929; #mol/liter
perc_glc=0.10;


l=cobra.io.read_sbml_model("../models/lactobacillus.sbml")
l.solver.configuration.timeout = 1 * 1000;
l.reactions.EX_lac_L_e_.id='EX_lac_e_';
l.name="Lactobacillus Brevis";
l.reactions.EX_o2_e_.lower_bound=-2.5;
l.reactions.EX_glc_e_.lower_bound=-18.5;
l.reactions.EX_nh4_e_.lower_bound=-1000;
ll=dmmm.Organism(l,vol*0.1,lambda x: x/30+1/(7*24))
ll.exchange_filter=lambda x: x.id!="added_biomass_sink";
ll.set_secondary_objective(l.reactions.EX_lac_e_.flux_expression,variance=0.02);
ll.michaelis_menten['EX_glc_e_']=0.1/MM['glc']*1000;
ll.michaelis_menten['EX_o2_e_']=0.005;
ll.inhibition['EX_glc_e_']={'EX_lac_e_':20000, 'EX_etoh_e_':50*220};

y=cobra.io.read_sbml_model("../models/yeast.sbml")
y.name="Sacharomyces Cerevisiæ";
y.reactions.EX_o2_e_.lower_bound=-2.5;
y.reactions.EX_glc_e_.lower_bound=-18.5;
yy=dmmm.Organism(y,500.0/200*vol,lambda x: x/27+1/(7*24)) # A cell dies: After 27 cell divisions or after two weeks of live
yy.inhibition['EX_glc_e_']={'EX_etoh_e_':220, 'EX_lac_e_':45}; #mmol/l
yy.set_secondary_objective(-y.reactions.EX_etoh_e_.flux_expression);
yy.michaelis_menten['EX_glc_e_']=0.5/MM['glc']*1000;
yy.michaelis_menten['EX_o2_e_']=0.005;

#co=[ll] #[ll,yy]
co=[yy] #[ll,yy]
#co=[ll,yy]

d=dmmm.DMMM(co)


d.medium.volume=lambda x: vol;

d.medium.metabolites["EX_h2o_e_"]=vol*mol_vol*1000;

#from paper
d.medium.metabolites["EX_glc_e_"]=vol*1000;
d.medium.metabolites["EX_nh3_e_"]=vol*(37.85+600);
d.medium.metabolites["EX_nh4_e_"]=vol*(37.85+600+4000);
d.medium.metabolites["EX_so4_e_"]=vol*(37.86 + 2.03 + 20 + 10 + 150);
d.medium.metabolites["EX_na1_e_"]=vol*(2.05 + 1.5);
d.medium.metabolites["EX_pi_e_"]=vol*(0.86 + 6.83 + 20+250);

d.medium.metabolites["EX_ribflv_e_"]=vol*(0.1); #200
d.medium.metabolites["EX_pnto_R_e_"]=vol*(1); #2000
d.medium.metabolites["EX_pydxn_e_"]=vol*(0.2); #400
d.medium.metabolites["EX_ribflv_e_"]=vol*(0.1); #200
d.medium.metabolites["EX_thm_e_"]=vol*(0.1); #200
d.medium.metabolites["EX_Nbfortyr_e_"]=vol*(0.1); #200

#From wikipedia
d.medium.metabolites["EX_btn_e_"]=vol*(0.1); #200

#From minrxns for yeast
for id in ['EX_zymst_e_', 'EX_ga6p_e_', 'EX_tre_e_']:
	d.medium.metabolites[id]=vol*100;

#minrxns for l.p.
for id in ['EX_2hxic_L_e_', 'EX_cys_L_e_', 'EX_ile_L_e_', 'EX_phe_L_e_', 'EX_pro_L_e_', 'EX_val_L_e_', 'EX_glu_L_e_', 'EX_pro_L_e_', 'EX_pydam_e_', 'EX_leu_e_']:
	d.medium.metabolites[id]=vol*30;

for id in ['EX_phe_L_e_', 'EX_leu_e_']:
	d.medium.metabolites[id]=vol*400;

#manual for l.p.
for id in ['EX_gln_L_e_', 'EX_ile_L_e_', 'EX_leu_L_e_', 'EX_met_L_e_', 'EX_nac_e_', 'EX_pydam_e_', 'EX_val_L_e_']:
	d.medium.metabolites[id]=vol*50;
d.medium.metabolites['EX_gln_L_e_']=vol*800;

d.medium.metabolites["EX_o2_e_"]=vol*0.5;

import matplotlib.pyplot as pp;

def plot_met(met):
	pp.plot(r.get_times(),r.get_metabolite("EX_"+met+"_e_"));

def plot_pop(idx):
	pp.plot(r.get_times(),r.get_population(idx));

def get_densities_at_point(idx):
	mets=d.split_state(r.state[idx])[1];
	densities=d.medium.list_to_densities(mets,vol);
	return densities;

def organism_at_point(org,idx):
	densities=get_densities_at_point(idx);
	o=d.organisms[org];
	o.adapt_boundaries(densities);
	return o;

def get_mass_densities(met, M):
	return [ r*M for r in r.get_metabolite(met) ];

def get_metabolite_density(r,met):
	return [val*MM[met]/1000/vol for val in r.get_metabolite("EX_"+met+"_e_")];

def generate_plots(r,basename,ext="pdf"):
	mets=["glc","etoh","co2","ac","lac","o2"];
	mets2=["nh3","so4","pi"];
	times=r.get_times();
	pp.clf();
	for met in mets:
		pp.semilogy(times,[val*MM[met]/1000/vol for val in r.get_metabolite("EX_"+met+"_e_")]);
	pp.legend(mets);
	pp.xlabel('time [h]');
	pp.ylabel('metabolite concentration [g/l]');
	pp.savefig(basename+"_metabolites."+ext);

	pp.clf();
	for i in range(0,len(r.dmmm.organisms)):
		pp.semilogy(times,[p/vol for p in r.get_population(i)]);
	pp.legend([o.model.name for o in r.dmmm.organisms]);
	pp.xlabel('time [h]');
	pp.ylabel('population [g/l]');
	pp.savefig(basename+"_populations."+ext);

import numpy as np;

def used_up_metabolites(r, threshold=1, fro=0, to=-1):
	s1=r.dmmm.split_state(r.state[fro])[1];
	s2=r.dmmm.split_state(r.state[to])[1];
	z1=s1<threshold;
	z2=s2<threshold;
	z=np.logical_and(z2,np.logical_not(z1));
	return [ r.dmmm.medium.lut[i] for i,v in enumerate(z) if v];


#r=d.simulate(25);

#import IPython; IPython.embed()
#exit()
rs=[]
concs=[0,22,44,66]
for conc in concs:
	d.medium.metabolites["EX_lac_e_"]=vol*conc;
	r=d.simulate(40);
	generate_plots(r, "yeast_adaption/"+str(conc)+"lac");
	rs.append(r);

pp.clf()
leg=[]
for r,c in zip(rs,concs):
	pp.plot(r.get_times(),get_metabolite_density(r,"glc"));
	pp.plot(r.get_times(),get_metabolite_density(r,"etoh"));
	leg.extend(["glc_"+str(c),"etoh_"+str(c)]);
pp.legend(leg);
pp.xlabel('time [h]');
pp.ylabel('metabolite densities [g/l]');
pp.savefig("yeast_adaption/similar_plot.pdf");

