import cobra;
import functools;
from collections import Counter;
import scipy.integrate;
import numpy as np;
import numbers;
import sys;

class DefaultDict(dict):
	def __init__(self, default):
		self.default=default;

	def __getitem__(self, k):
		return self.get(k);

	def get(self, key, default=None):
		if not default:
			default=self.default
		return super(DefaultDict,self).get(key,default);

class Organism:
	def __init__(self, model, population, mortality):
		self.mortality=mortality;
		self.population=population;
		self.model=model;
		self.inhibition={};
		self.michaelis_menten=DefaultDict(0.1);
		self.exchange_filter=None;
		self.initial_lb={r.id:r.lower_bound for r in self.exchanges()};
		self.secondary_objective=None;

	def density_to_lower_bound(self, densities, id):
		s=densities[id];
		vm=self.initial_lb[id];	#This is usually negative
		k=self.michaelis_menten[id];
		tox=1.0;
		for oid,inhibition in self.inhibition.get(id,{}).items():
			tox=tox*1/(1+densities[oid]/inhibition);
		return vm*s/(k+s)*tox;		#Therefore, this is usually too

	def lower_bounds(self, densities):
		return {k:self.density_to_lower_bound(densities, k) for k in densities.keys()};

	def adapt_boundaries(self, densities):
		for id,r in [(r.id,r) for r in self.exchanges()]:
			r.lower_bound=self.density_to_lower_bound(densities,id);

	def exchanges(self):
		if self.exchange_filter:
			return filter(self.exchange_filter,self.model.exchanges);
		return self.model.exchanges;

	def set_secondary_objective(self,expr,direction='max'):
		self.condition=self.model.problem.Constraint(self.model.objective.expression);
		self.secondary_objective=self.model.problem.Objective(expr, direction=direction);

	def turn_off_addcon(self):
		self.model.remove_cons_vars(self.condition);

	def set_addcon(self,val):
		if not self.condition.ub or val > self.condition.lb:
			self.condition.ub=val;
			self.condition.lb=val;
		else:
			self.condition.lb=val;
			self.condition.ub=val;
		self.model.add_cons_vars(self.condition);

	def get_mortality(self,sol):
		if isinstance(self.mortality,numbers.Number):
			return self.mortality
		else:
			return self.mortality(sol);

	def fba(self,medium,densities):
		self.adapt_boundaries(densities);
		sol=None;
		if self.secondary_objective:
			val=self.model.slim_optimize();
			if self.model.solver.status != 'optimal':
				return (-self.get_mortality(self.model.optimize()),[0]);
			old_obj=self.model.objective;
			self.set_addcon(val);
			self.model.objective=self.secondary_objective;
			sol=self.model.optimize();
			self.model.objective=old_obj;
			self.turn_off_addcon()
			sol.objective_value=val;
		else:
			sol=self.model.optimize();

		if sol.status!='optimal':
			return (-self.get_mortality(sol),[0]);

		return (
			sol.objective_value-self.get_mortality(sol),
			medium.dict_to_list(sol[[r.id for r in self.exchanges()]])
		);

class Medium:

	def __init__(self, organisms):
		ids=map(lambda x: {r.id for r in x.exchanges()},organisms);
		medium_ids=functools.reduce(lambda x,y: x.union(y), ids);
		self.metabolites=dict.fromkeys(medium_ids, 0.0);
		self.lut=[id for id in self.metabolites.keys()];
		self.volume=lambda x: sum(x);

	def get_volume(self):
		return self.volume(self.metabolites);

	def densities(self):
		volume=self.get_volume();
		return {k:v/total for k,v in self.metabolites.items()};

	def list_to_dict(self, l):
		return {k:v for k,v in zip(self.lut,l)};

	def list_to_densities(self, l, total):
		return {k:v/total for k,v in zip(self.lut,l)};

	def dict_to_list(self, d):
		return [d.get(k,0) for k in self.lut];

class Result:

	def __init__(self,dmmm):
		self.dmmm=dmmm;
		self.state=[];
		self.times=[];

	def append(self,t,y):
		y=y.copy()
		y[y<0]=0;
		self.times.append(t);
		self.state.append(y);

	def get_population(self, idx):
		return [self.dmmm.split_state(s)[0][idx] for s in self.state];

	def get_metabolite(self, name):
		idx=self.dmmm.medium.lut.index(name);
		return [self.dmmm.split_state(s)[1][idx] for s in self.state];

	def get_times(self):
		return self.times;

	def get_metabolites_at_step(self, step):
		return self.dmmm.medium.list_to_dict(self.dmmm.split_state(self.state[step])[1]);

class DMMM:
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
		neg=y<0;
		self.buf[:]=[max(0,v) for v in y];
		(populations,metabolites)=self.split_state(self.buf);
		(growth,fluxes)=self.split_state(self.ret);
		self.ret[:]=0;
		
		total=self.medium.get_volume();
		self.densities.update(self.medium.list_to_densities(metabolites,total));
		
		for idx,o in enumerate(self.organisms):
			if populations[idx]>0: #No need to calculate the fba on dead populations
				(bio,met)=o.fba(self.medium, self.densities);
				growth[idx]=bio*populations[idx];
				fluxes[:]=fluxes+[m*populations[idx] for m in met];
		self.ret[neg]=[max(0,i) for i in self.ret[neg]];
		return self.ret;

	def get_ode(self):
		ode=scipy.integrate.ode(self.func);
		l=[];
		l.extend([v.population for v in self.organisms]);
		l.extend([self.medium.metabolites[id] for id in self.medium.lut]);
		ode.set_initial_value(l,0);
		#ode.set_f_params(self);
		ode.set_integrator('dopri5',atol=1e-5,rtol=1e-5);
		return ode;

	def simulate(self, t):
		ode=self.get_ode();
		r=Result(self);
		ode.set_solout(lambda t,y: r.append(t,y));
		ode.integrate(t);
		return r;

