import cobra;
import functools;
from collections import Counter;
import scipy.integrate;
import numpy as np;

class Organism:
	mortality=0.0;
	population=1.0;
	alpha=1.0;
	model=None;
	exchange_filter=None;
	initial_lb={};

	def __init__(self, model, population, mortality):
		self.mortality=mortality;
		self.population=population;
		self.model=model;
		self.initial_lb={r.id:r.lower_bound for r in self.exchanges()};
		retval=np.zeros(len(model.metabolites));

	def density_to_lower_bound(self, densities, id):
		s=densities[id];
		vm=self.initial_lb[id];	#This is usually negative
		k=0.1;
		return vm*s/(k+s);		#Therefore, this is usually too

	def adapt_boundaries(self, densities):
		for id,r in [(r.id,r) for r in self.exchanges()]:
			r.lower_bound=self.density_to_lower_bound(densities,id);
	
	def exchanges(self):
		if self.exchange_filter:
			return filter(self.exchange_filter,self.model.exchanges);
		return self.model.exchanges;

	def fba(self,medium,densities):
		self.adapt_boundaries(densities);
		sol=self.model.optimize();
		if sol.status!='optimal':
			return (0,0);
		assert(sol.objective_value >= 0);

		return (
			sol.objective_value*self.alpha-self.mortality,
			medium.dict_to_list(sol[[r.id for r in self.exchanges()]])
		);

class Medium:
	metabolites={};
	lut=[];

	def __init__(self, organisms):
		ids=map(lambda x: {r.id for r in x.exchanges()},organisms);
		medium_ids=functools.reduce(lambda x,y: x.union(y), ids);
		self.metabolites=dict.fromkeys(medium_ids, 0.0);
		self.lut=[id for id in self.metabolites.keys()];

	def densities(self):
		total=sum(self.metabolites);
		return {k:v/total for k,v in self.metabolites.items()};
		
	def list_to_dict(self, l):
		return {k:v for k,v in zip(self.lut,l)};
		
	def list_to_densities(self, l, total):
		return {k:v/total for k,v in zip(self.lut,l)};
		
	def dict_to_list(self, d):
		return [d.get(k,0) for k in self.lut];

class Result:
	dmmm=None;
	state=[];
	times=[];
	
	def __init__(self,dmmm):
		self.dmmm=dmmm;
	
	def append(self,t,y):
		self.times.append(t);
		self.state.append(y);
	
	def get_population(self, idx):
		idx=self.dmmm.medium.lut.index(name);
		return {t:self.dmmm.split_state(self.state[i])[1][idx] for i,t in enumerate(self.times)};
	
	def get_metabolite(self, name):
		idx=self.dmmm.medium.lut.index(name);
		return {t:self.dmmm.split_state(self.state[i])[1][idx] for i,t in enumerate(self.times)};
	
	def get_metabolites_at_step(self, step):
		return self.dmmm.medium.list_to_dict(self.dmmm.split_state(self.state[step])[1]);

class DMMM:
	organisms=[];
	medium=None;
	ret=None;
	buf=None;
	densities=None;

	def __init__(self, organisms):
		self.organisms=organisms;
		self.medium=Medium(organisms);
		self.ret=np.zeros(self.state_length());
		self.buf=np.zeros(self.state_length());
		self.densities={k:0 for k in self.medium.lut};

	def state_length(self):
		return len(self.medium.metabolites)+len(self.organisms);

	def simulate(self, delta_t):
		self.adapt_boundaries();
		for organism in self.organisms:
			dM=organism.simulate(delta_t);
			for k,v in dM.items():
				self.medium.metabolites[k]+=v;
	
	def split_state(self,s):
		return (s[0:len(self.organisms)],s[len(self.organisms):]);
	
	def func(self, t, y):
		self.buf[:]=[max(0,v) for v in y];
		(populations,metabolites)=self.split_state(self.buf);
		(growth,fluxes)=self.split_state(self.ret);
		self.ret[:]=0;
		
		#print("State");
		#print(populations);
		#print(self.medium.list_to_dict(metabolites));
		total=sum(metabolites);
		self.densities.update(self.medium.list_to_densities(metabolites,total));
		
		for idx,o in enumerate(self.organisms):
			if populations[idx]>0: #No need to calculate the fba on dead populations
				(bio,met)=o.fba(self.medium, self.densities);
				growth[idx]=bio*populations[idx];
				fluxes[:]=fluxes*populations[idx]+met;
		#print("Derivative");
		#print(growth)
		#print(self.medium.list_to_dict(fluxes));
		#assert(False);
		return self.ret;

	def get_ode(self):
		ode=scipy.integrate.ode(self.func);
		l=[];
		l.extend([v.population for v in self.organisms]);
		l.extend([self.medium.metabolites[id] for id in self.medium.lut]);
		ode.set_initial_value(l,0);
		#ode.set_f_params(self);
		ode.set_integrator('dopri5')
		return ode;

	def simulate(self, t):
		ode=self.get_ode();
		r=Result(self);
		ode.set_solout(lambda t,y: r.append(t,y.copy()));
		ode.integrate(t);
		return r;

