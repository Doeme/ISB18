import cobra;
import functools;

class Organism:
	mortality=0.0;
	population=1.0;
	alpha=1.0;
	model=None;
	def __init__(self, model, population, mortality):
		self.mortality=mortality;
		self.population=population;
		self.model=model;

	def density_to_lower_bound(density):
		return density*0.25;

	def adapt_boundaries(self, medium):
		for met,d in medium.densities().items():
			model.reactions.get_by_id(met).lower_bound=density_to_lower_bound(d);

	def simulate(self,delta_t):
		sol=self.model.optimize();
		assert(sol.status=='optimal');
		metabolic_fluxes=dict(sol[[r.id for r in self.model.exchanges]]);
		mf=metabolic_fluxes;
		biomass=sol.objective_value;

		growth_factor=(biomass/self.alpha-self.mortality);
		if abs(growth_factor) < sys.float_info.epsilon: #Prevent division by zero
			new_population=self.population;
			met_factor = self.population*delta_t;
		else:
			new_population=self.population*exp(growth_factor*delta_t);
			met_factor = self.population/growth_factor * (exp(growth_factor*delta_t)-1);

		self.population=new_population;
		return {k: v*met_factor for k,v in mf.items()};

class Medium:
	metabolites={};

	def __init__(self, organisms):
		models=map(lambda x: x.model,organisms);
		ids=map(lambda x: {r.id for r in x.exchanges},models);
		medium_ids=functools.reduce(lambda x,y: x.union(y), ids);
		self.metabolites=dict.fromkeys(medium_ids, 0.0);

	def densities(self):
		total=sum(self.medium);
		return {k:v/total for k,v in self.medium.items()};

class DMMM:
	organisms=[];
	medium=None;

	def __init__(self, organisms):
		self.organisms=organisms;
		self.medium=Medium(organisms);

	def simulate(self, delta_t):
		for organism in self.organisms:
			organism.adapt_boundaries(self.medium);
			dM=organism.simulate(delta_t);
			for k,v in dM.items():
				medium.metabolites[k]+=v;
