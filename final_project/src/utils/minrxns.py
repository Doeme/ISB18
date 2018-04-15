import cobra

def minrxns(m, fac=0.1, lb=True, ub=False, reactions=None):
	if not reactions:
		reactions=m.exchanges;
	with m:
		sol=m.slim_optimize();
		bio=None;
		for rxn in m.objective.variables:
			if rxn.name in m.reactions:
				bio=rxn.name;
				break;
		m.reactions.get_by_id(bio).lower_bound=sol*fac;
	
		expr=None;
		rxns=[];
		for rxn in reactions:
			n=rxn.id;
			var=m.problem.Variable(n+'_knockout',type='binary');
			lbv=rxn.lower_bound;
			ubv=rxn.upper_bound;
			if not expr:
				expr=var;
			else:
				expr=expr+var;
			m.add_cons_vars(var);
			if ub:
				cu=m.problem.Constraint(rxn.flux_expression-ubv*var,ub=0);
				m.add_cons_vars(cu);
			if lb:
				cl=m.problem.Constraint(rxn.flux_expression-lbv*var,lb=0);
				m.add_cons_vars(cl);
			rxns.append((rxn,var))
	
		m.objective=m.problem.Objective(expr, direction='min');
		sol=m.optimize();
		important=[r[0].id for r in rxns if r[1].primal > 0.5];
		unimportant=[r[0].id for r in rxns if r[1].primal < 0.5];
	return (important,unimportant);
