import cobra;
import cobra.test;
m = cobra.test.create_test_model("textbook")
#m=cobra.io.load_json_model('../ecolicore.json');
bio=m.reactions.get_by_id('Biomass_Ecoli_core');
succ=m.reactions.get_by_id('EX_succ_e');
optbio=m.slim_optimize();
m.objective=succ;
optsucc=m.slim_optimize();

import cobra.flux_analysis.variability as fvam
fva=fvam.flux_variability_analysis

with m:
	bio.lower_bound=0.9*optbio;
	bio.upper_bound=0.9*optbio;
	fvai=fva(m);

with m:
	succ.lower_bound=0.9*optsucc;
	succ.upper_bound=0.9*optsucc;
	fvaii=fva(m);

import sys;
eps=sys.float_info.epsilon*100;

def interesting(id):
	r=m.reactions.get_by_id(id);
	if r.boundary or len(r.compartments)>1:
		return False;
	if abs(fvai['minimum'][id])<eps and abs(fvai['maximum'][id])<eps and abs(fvaii['minimum'][id])<eps and abs(fvaii['maximum'][id])<eps:
		return False;
	return True;

fvai_i=fvai.select(interesting);
fvaii_i=fvaii.select(interesting);

f=fvai_i.join(fvaii_i,lsuffix="_i",rsuffix="_ii");



def non_overlap(a,b,c,d):
	return (d>a-eps or c<b+eps);

def non_overlap_id(fva,id):
	return non_overlap(fva['maximum_i'][id],fva['minimum_i'][id],fva['maximum_ii'][id],fva['minimum_ii'][id]);

ff=f.select(lambda x: non_overlap_id(f,x));
print(", ".join([i for i in ff.index]));

mini=min(ff.min())
maxi=max(ff.max())
print("min/max: {}/{}".format(mini,maxi));

h=open("barplot.tex", "w+");

for n,a,b,c,d in f.itertuples():
	a=min(a,maxi);
	c=min(c,maxi);
	b=max(b,mini);
	d=max(d,mini);
	col='black';
	if(non_overlap(a,b,c,d)):
		col='red';
	
	print("{{\\color{{{}}}\\verb|{}|}} & \\barplot".format(col,n), file=h, end="",sep="")
	for tupl in [(a,b), (c,d)]:
		if(abs(tupl[1]-tupl[0])<1):
			print("{{({0:.4f},0) circle(0.25ex)}}".format((tupl[0]+tupl[1])/2),file=h,end="",sep="")
		else:
			print("{{({0:.4f},0) -- ({1:.4f},0)}}".format(tupl[0],tupl[1]),file=h,end="",sep="")
	
	print("\\\\",file=h);

h.close()
