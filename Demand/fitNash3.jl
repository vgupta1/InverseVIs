## File contains the fitting funcs with a parametric expansion...

using JuMP
using Gurobi
using Distributions
using Roots

#represent the desired funciton as a combination of features
#featureFuncs is itself a collection of function pointers
# phi : (prices, gdp) -> Real

#Cheapest way possible
function genFeaturesPoly(deg)
	features = Array(Any, binomial(deg + 3, deg))
	indx = 1
	for p1deg = 0:deg
		for p2deg = 0:deg
			for gdpdeg = 0:deg
				if p1deg + p2deg + gdpdeg <= deg
					features[indx] = (p1, p2, gdp) -> p1^p1deg * p2^p2deg * gdp^gdpdeg
					indx+=1
				end
			end
		end
	end
	features
end

function genFeaturesExp( rng )
	features = [(p1, p2, gdp)-> exp(r * p1) for r in rng]
	append!(features, [(p1, p2, gdp) -> exp(r * p2) for r in rng])
	append!(features, [(p1, p2, gdp) -> exp(r * gdp) for r in rng])
	append!(features, [(p1, p2, gdp) -> p1, (p1, p2, gdp) -> p2, (p1, p2, gdp) -> gdp, (p1, p2, gdp)-> 1])
	features
end

createfun(flambdas, features) = 
	(prices, gdp) -> sum([flambdas[i] * features[i](prices..., gdp) for i=1:length(flambdas) ])

function setUpFitting(features)
	m = Model(solver=GurobiSolver(OutputFlag=false))
#	m = Model()
	N=length(features)
	@defVar(m, vlambdas1[1:N])
	@defVar(m, vlambdas2[1:N])
	return m, vlambdas1, vlambdas2
end

#Gens the exp features and coefficients imultaneously to constrain increasingness etc.
function setUpFitting( rng::Range{Float64} )
	N = 3*length(rng) + 4 # add linears and the constants
	features = genFeaturesExp(rng)
	m = Model(solver=GurobiSolver(OutputFlag=false))
	@defVar(m, vlambdas1[1:N])
	@defVar(m, vlambdas2[1:N])

	#require that MR(infty) -> -infty
	# for (ix, r) = enumerate(rng)
	# 	if r >= 0
	# 		@addConstraint(m, vlambdas1[ix] <= 0)
	# 		@addConstraint(m, vlambdas2[ix + length(rng)] <= 0)
	# 	end
	# 	# if r >= 0  #increasing functions
	# 	# 	@addConstraint(m, vlambdas1[ix] == 0)  #MR1 should be decreasing in p1
	# 	# 	@addConstraint(m, vlambdas1[ix + length(rng)] >=0 )  #MR1 should be increasing in p2
	# 	# 	@addConstraint(m, vlambdas2[ix] >= 0)  #MR2 should be increasing in p1
	# 	# 	@addConstraint(m, vlambdas2[ix + length(rng)] == 0)  #MR2 should be decreasing in p2
	# 	# elseif r <= 0  #decreasing functions
	# 	# 	@addConstraint(m, vlambdas1[ix] >= 0)   #MR1 should be decreasing in p1
	# 	# 	@addConstraint(m, vlambdas1[ix + length(rng)] == 0)   #MR1 should be increasing in p2
	# 	# 	@addConstraint(m, vlambdas2[ix] == 0)  #MR2 should be increasing in p1
	# 	# 	@addConstraint(m, vlambdas2[ix + length(rng)] >=0)  #MR2 should be decreasing in p2
	# 	# end
	# end

	#same increasing/decreasing for the linear terms
	# @addConstraint(m, vlambdas1[3 *length(rng) + 1] <= 0)
	# @addConstraint(m, vlambdas2[3 *length(rng) + 2] <= 0)
	# @addConstraint(m, vlambdas1[3 *length(rng) + 2] >= 0)
	# @addConstraint(m, vlambdas2[3 *length(rng) + 1] >= 0)

	#disallow the weirdness in the functional form
	for ix = 1:length(rng)
		@addConstraint(m, vlambdas1[ix + 2 *length(rng)] == 0)
		@addConstraint(m, vlambdas2[ix + 2 *length(rng)] == 0)
	end
	return m, vlambdas1, vlambdas2, features
end

#insist that the fitted function is decreasing on the data points for mean values of other stuff		 
function addDecreasingCnsts(m, vlambdas, features, otherP, meanGDP, price_data, pmax, ID)
	price_data = sort(price_data)
	sort_prices(myP, yourP) = ID == :1 ? (myP, yourP) : (yourP, myP)
	for j = 1:(length(price_data)- 1)
		prices_j = sort_prices(price_data[j], otherP)
		prices_jp = sort_prices(price_data[j+1], otherP)
		@addConstraint(m, sum{features[i](prices_j..., meanGDP) * vlambdas[i], i=1:length(vlambdas)} >=
		 					sum{features[i](prices_jp..., meanGDP) * vlambdas[i], i=1:length(vlambdas)} )
	end

	#add one last one for pmax in case its not achieved...
	prices_j = sort_prices(price_data[length(price_data)], otherP)
	prices_jp = sort_prices(pmax, otherP)
	@addConstraint(m, sum{features[i](prices_j..., meanGDP) * vlambdas[i], i=1:length(vlambdas)} >=
	 					sum{features[i](prices_jp..., meanGDP) * vlambdas[i], i=1:length(vlambdas)} )

end

#insist that the fitted function is positive at zero for mean values
function addPositivity(m, vlambdas, features, otherP, myP, meanGDP, ID; TOL=0)
	if ID == :1
		@addConstraint(m, sum{features[i](myP, otherP, meanGDP) * vlambdas[i], i=1:length(vlambdas)} >= TOL)
	else
		@addConstraint(m, sum{features[i](otherP, myP, meanGDP) * vlambdas[i], i=1:length(vlambdas)} >= TOL)
	end
end


#insist that the MR at the mean point should be equal to given  quantity for each
function normalize(m, vlambdas1, vlambdas2, features, meanPs, meanGDP, MR1, MR2)
	N = length(vlambdas1)
	@addConstraint(m, MR1 == sum{features[i](meanPs..., meanGDP)*vlambdas1[i], i=1:N})
	@addConstraint(m, MR2 == sum{features[i](meanPs..., meanGDP)*vlambdas2[i], i=1:N})
end

function normalize(m, vlambdas, features, prices, error, MR)
	N = length(vlambdas)
	@addConstraint(m, MR == sum{features[i](prices..., error) * vlambdas[i], i=1:N})
end



#adds the constraints to the model for the residual for this observation
function addResids(m, vlambdas1, vlambdas2, prices_t, GDP_t, pmax, features)
	#add some supplementary variables for ease:
	x = [prices_t GDP_t]; N = length(vlambdas1)
	@defVar(m, MR[1:2])
	@addConstraint(m, MR[1] == sum{features[i](prices_t..., GDP_t)*vlambdas1[i], i=1:N})
	@addConstraint(m, MR[2] == sum{features[i](prices_t..., GDP_t)*vlambdas2[i], i=1:N})

	@defVar(m, resid >= 0)
	@defVar(m, y[1:2] >= 0)
	@addConstraint(m, y[1] >= MR[1])
	@addConstraint(m, y[2] >= MR[2])
	@addConstraint(m, pmax * (y[1] + y[2]) - prices_t[1] * MR[1] - prices_t[2] * MR[2] <= resid)
	return resid
end

#calc an out of sample residual.  may be a better way to do this.
function calcResid(flambdas1, flambdas2, prices, GDP, pmax, features)
	N = length(flambdas1); MR = Array(Float64, 2)
	MR[1] = sum([flambdas1[i]*features[i](prices..., GDP) for i = 1:N])
	MR[2] = sum([flambdas2[i]*features[i](prices..., GDP) for i = 1:N])

	y = Array(Float64, 2)
	y[1] = max(0, MR[1])
	y[2] = max(0, MR[2])

	return pmax * (y[1] + y[2]) - prices[1] * MR[1] - prices[2] * MR[2]
end

function bestResponsePrice(flambdas, features, prices, error, pmax, ID)
	fMR = createfun(flambdas, features)
	if ID == :1
		# pstar = newton(p-> fMR([p, prices[2]], error), prices[1])
		pstar = fzero(p-> fMR([p, prices[2]], error), prices[1])
		#pstar = customRoot(p-> fMR([p, prices[2]], error), pmax)  #works really well with exp features
	else	
		# pstar = newton(p-> fMR([prices[1], p], error), prices[2])
		pstar = fzero(p-> fMR([prices[1], p], error), prices[2])
		# pstar = customRoot(p-> fMR([prices[1], p], error), pmax)   #works really well with exp features
	end
	# assert( pstar >= 0 )
	pstar = max(0, pstar)
	return min(pstar, pmax)
end

#assumes f(0) > 0 and we have decreasingness
function customRoot(f, pmax; MAX_BRACKET = 100)
	assert( f(1e-3) > 0 )  # this should be because of the positivity issue
	if f(pmax) > 0
		return pmax
	end
	return fzero(f, [1e-3, pmax]) 
end

#solve for equilibrium prices.  True model is log with a random, firm dependent effect
function solveNashPrices(flambdas1, flambdas2, features, error, pmax; MAX_ITERS = 20, TOL=1e-6, trace=false, 
		prices=fill(convert(Float64, pmax / 2), 2))
	dist = 0
	for i = 1:MAX_ITERS
		prices_old = copy(prices)
		prices[1] = bestResponsePrice(flambdas1, features, prices_old, error, pmax, :1)
		prices[2] = bestResponsePrice(flambdas2, features, prices_old, error, pmax, :2)
		dist = norm(prices - prices_old)
		if trace
			println([prices; dist]')
		end
		if dist < TOL
			return prices
		end
	end
	println("Warning: Max Iters Reached \t $dist")
	return prices
end


