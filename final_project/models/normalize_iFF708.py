import cobra;

#SLFxtO SLFxtI
m=cobra.io.read_sbml_model("iFF708.xml");

lut={'ETH':'etoh','SLF':'so4','NA':'na1','RFLAV':'ribflv'};

for r in m.exchanges:
	n=r.id;
	if n[-1]=='O':
		on=r.id[:-1]+'I';
		orx=m.reactions.get_by_id(on);
		r.lower_bound=r.lower_bound-orx.upper_bound;
		r.upper_bound=r.upper_bound-orx.lower_bound;
		orx.remove_from_model();
		met=n[:-3];
		met=lut.get(met,met);
		r.id="EX_"+met.lower()+"_e_";

cobra.io.write_sbml_model(m,"iFF708_fixed.xml");
