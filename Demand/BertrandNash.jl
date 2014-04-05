# Bertrand Nash example
# Generates random equilibrium data
function bestResponsePrice(thetas, intercept, prices, pmax, FIRM_ME)
	if FIRM_ME == :1
		me = 1; other = 2
	else
		me = 2; other = 1
	end
	p = (intercept + thetas[other] * prices[other]) / 2 / abs(thetas[me])
	if p <= 0
		return 0
	elseif p >= pmax
		return pmax
	else
		return p
	end
end

#solve for equilibrium prices.  True model is linear with random effects
function solveNashPrices(thetas1, thetas2, intercept1, intercept2, pmax; MAX_ITERS = 20, TOL=1e-6, trace=false)
	prices = fill(convert(Float64, pmax), 2)
	for i = 1:MAX_ITERS
		prices_old = copy(prices)
		prices[1] = bestResponsePrice(thetas1, intercept1, prices, pmax, :1)
		prices[2] = bestResponsePrice(thetas2, intercept2, prices, pmax, :2)
		dist = norm(prices - prices_old)
		if trace
			println( prices', dist )
		end
		if dist < TOL
			return prices
		end
	end
	error("Max Iters reached")
end

solveNashPricesFit(thetas1, thetas2, intercept1, intercept2, GDP, pmax) = 
	solveNashPrices(thetas1[1:2], thetas2[1:2], thetas1[3] * GDP + intercept1, thetas2[3] * GDP + intercept2, pmax)

function genErrors()
	eps1, eps2, eps3 = rand(Normal(5, .2), 3)
	eps3 += eps1 + eps2
	return eps1, eps2, eps3
end

