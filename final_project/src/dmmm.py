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
	initial_lb={};
	def __init__(self, model, population, mortality):
		self.mortality=mortality;
		self.population=population;
		self.model=model;
		self.initial_lb={r.id:r.lower_bound for r in model.exchanges};

	def density_to_lower_bound(self, densities, id):
		s=densities[id];
		vm=abs(self.initial_lb[id]);
		k=1;
		return -vm*s/(k+s);

	def adapt_boundaries(self, densities):
		for id,r in [(r.id,r) for r in self.model.exchanges]:
			r.lower_bound=self.density_to_lower_bound(densities,id);

	def fba(self,densities):
		self.adapt_boundaries(densities);
		sol=self.model.optimize();
		if sol.status!='optimal':
			return (0,{});
		assert(sol.objective_value >= 0);

		return (
			sol.objective_value*self.alpha-self.mortality,
			dict(sol[[r.id for r in self.model.exchanges]])
		);

class Medium:
	metabolites={};

	def __init__(self, organisms):
		models=map(lambda x: x.model,organisms);
		ids=map(lambda x: {r.id for r in x.exchanges},models);
		medium_ids=functools.reduce(lambda x,y: x.union(y), ids);
		self.metabolites=dict.fromkeys(medium_ids, 0.0);

	def densities(self):
		total=sum(self.metabolites);
		return {k:v/total for k,v in self.metabolites.items()};

class DMMM:
	organisms=[];
	medium=None;

	def __init__(self, organisms):
		self.organisms=organisms;
		self.medium=Medium(organisms);

	def adapt_boundaries(self):
		for organism in self.organisms:
			organism.adapt_boundaries(self.medium);

	def simulate(self, delta_t):
		self.adapt_boundaries();
		for organism in self.organisms:
			dM=organism.simulate(delta_t);
			for k,v in dM.items():
				self.medium.metabolites[k]+=v;

	def func(self, t, y):
		c=Counter({ k:0 for k in self.medium.metabolites.keys()});
		populations=y[0:len(self.organisms)];
		metabolites=y[len(self.organisms):];
		total=sum(metabolites);
		densities={k:metabolites[idx]/total for idx,k in enumerate(self.medium.metabolites.keys())};
		ret=np.zeros(len(y));
		for idx,o in enumerate(self.organisms):
			(bio,met)=o.fba(densities);
			ret[idx]=bio;
			c=c+Counter({k:v*y[idx] for k,v in met.items()});
		ret[len(self.organisms):]=[c[id] for id in self.medium.metabolites.keys()];
		return ret;

	def get_ode(self):
		ode=scipy.integrate.ode(self.func);
		l=[];
		l.extend([v.population for v in self.organisms]);
		l.extend([v for v in self.medium.metabolites.values()]);
		ode.set_initial_value(l,0);
		#ode.set_f_params(self);
		ode.set_integrator('dopri5')
		return ode;

	def simulate(self, t):
		ode=self.get_ode();
		l=[];
		ode.set_solout(lambda t,y: print(t,y)); #l.append((t,y)));
		ode.integrate(t);
		return l;
