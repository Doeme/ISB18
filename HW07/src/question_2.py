import cobra;

#m=cobra.io.load_matlab_model("../Model_iJO1366.mat")
m=cobra.io.load_json_model("../ecolicore.json")
asp_c=cobra.Metabolite('asp_c')
ala_c=cobra.Metabolite('ala_c')
oxo_c=cobra.Metabolite('oxo_c')
hpa_c=cobra.Metabolite('hpa_c')
glu_c=cobra.Metabolite('glu_c')

aat=cobra.Reaction('AAT')
pand=cobra.Reaction('PAND')
gabt=cobra.Reaction('GABT')
hpdh=cobra.Reaction('HPDH')
hpa_ex=cobra.Reaction('HPA_Ex')

akg_c=m.metabolites.get_by_id('akg_c');
oaa_c=m.metabolites.get_by_id('oaa_c');
co2_c=m.metabolites.get_by_id('co2_c');
nadph_c=m.metabolites.get_by_id('nadph_c');
nadp_c=m.metabolites.get_by_id('nadp_c');
aat.add_metabolites({oaa_c:-1,asp_c:1,glu_c:-1,akg_c:1})
pand.add_metabolites({asp_c:-1,ala_c:1,co2_c:1});
gabt.add_metabolites({ala_c:-1,akg_c:-1,glu_c:1,oxo_c:1})
hpdh.add_metabolites({oxo_c:-1,nadph_c:-1,hpa_c:1,nadp_c:1})
hpa_ex.add_metabolites({hpa_c:-1});
hpa_ex.lower_bound=0;

m.add_metabolites([asp_c,ala_c,oxo_c,hpa_c,glu_c]);
m.add_reactions([aat,pand,gabt,hpdh,hpa_ex]);

oxex=m.reactions.get_by_id('EX_o2_e');

def sol_intakes(sol):
	exch=sol.fluxes[[ r.id for r in m.exchanges]];
	return exch.select(lambda x: exch[x]<0);

def sol_secretion(sol):
	exch=sol.fluxes[[ r.id for r in m.exchanges]];
	return exch.select(lambda x: exch[x]>0);

def sol_yield(sol):
	glc=sol.fluxes['EX_glc__D_e'];
	return abs(sol.objective_value/glc);

oxex.lower_bound=-1000;
sol=m.optimize();
print("i) Aerobic conditions yield (ex_hpa/ex_glc): {}".format(sol_yield(sol)))
oxex.lower_bound=0;
sol=m.optimize();
print("ii) Anaerobic conditions yield (ex_hpa/ex_glc): {}".format(sol_yield(sol)))

import numpy as np;
import itertools;

glcex=m.reactions.get_by_id('EX_glc__D_e');
oxr=np.linspace(-10,0,11);
glcr=np.linspace(-10,0,11);
points=itertools.product(oxr,glcr);
res=[]
for ox,glc in points:
	oxex.lower_bound=ox;
	oxex.upper_bound=ox;
	glcex.lower_bound=glc;
	glcex.upper_bound=glc;
	sol=m.optimize('maximize');
	res.append([ox,glc,sol.objective_value]);

import matplotlib.pyplot as plt
plt.scatter([p[1] for p in res],[p[2] for p in res],c=[p[0] for p in res])
#plt.fill_between();

ax=plt.axes();
ax.set_xlabel('Glucose Intake');
ax.set_ylabel('HPA Secretion');
cbar=plt.colorbar();
cbar.ax.set_ylabel('Oxygen Intake');

plt.grid();

plt.savefig("phaseplane.pdf");
