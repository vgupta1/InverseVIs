## File contains the fitting funcs


using JuMP
using Gurobi
using Distributions

include("BertrandNash.jl")  #file generates data

function setUpFitting()
	m = Model(solver=GurobiSolver(OutputFlag=false))
	@defVar(m, vthetas1[1:3])
	@defVar(m, vthetas2[1:3])
	# @addConstraint(m, vthetas1[2])
	# @addConstraint(m, vthetas2[1])
	@defVar(m, vintercept1 >= 0)
	@defVar(m, vintercept2 >= 0)
	#normalization choice
	@addConstraint(m, vthetas1[1] <= -1.)
	@addConstraint(m, vthetas2[2] <= -1.)
	return m, vthetas1, vthetas2, vintercept1, vintercept2
end

#adds the constraints to the modl for the residual for this observation
function addResids(m, thetas1, thetas2, intercept1, intercept2, prices_t, GDP_t, pmax)
	@defVar(m, y[1:2] >= 0)
	@defVar(m, resid >= 0)
	@defVar(m, d[1:2] >=0)  #doesn't seem to work if we insist these are positive...
	@addConstraint(m, d[1] == prices_t[1] * thetas1[1] + prices_t[2] * thetas1[2] + GDP_t * thetas1[3] + intercept1 )
	@addConstraint(m, d[2] == prices_t[1] * thetas2[1] + prices_t[2] * thetas2[2] + GDP_t * thetas2[3] + intercept2 )
	@addConstraint(m, y[1] >= d[1] + prices_t[1] * thetas1[1])
	@addConstraint(m, y[2] >= d[2] + prices_t[2] * thetas2[2])
	@addConstraint(m, pmax * (y[1] + y[2]) - prices_t[1] * d[1] - prices_t[2] * d[2] - 
					prices_t[1]^2 * thetas1[1] - prices_t[2]^2 * thetas2[2] <= resid )
	return resid
end

#calc an out of sample residual.  may be a better way to do this.
function calcResid(fitTheta1, fitTheta2, fitIntercept1, fitIntercept2, prices, GDP, pmax)
	m, vthetas1, vthetas2, vintercept1, vintercept2 = setUpFitting()
	for i = 1:3
		@addConstraint(m, vthetas1[i] == fitTheta1[i])
		@addConstraint(m, vthetas2[i] == fitTheta2[i])
	end
	@addConstraint(m, vintercept1 == fitIntercept1)
	@addConstraint(m, vintercept2 == fitIntercept2)
	resid = addResids(m, vthetas1, vthetas2, vintercept1, vintercept2, prices, GDP, pmax)
	@setObjective(m, Min, resid)
	status = solve(m); 
	if status != :Optimal
		println(prices, status)
		error(status)
	end

	return getObjectiveValue(m)
end

