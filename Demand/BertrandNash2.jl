# Bertrand Nash example
# Generates random equilibrium data
using Distributions

indxID(ID) = (ID == :1) ? (1, 2) : (2, 1)

function bestResponsePrice(thetas, prices, error, pmax, ID)
	me, you = indxID(ID)
	b = thetas[you]*prices[you] + thetas[3]*error + thetas[4]
	p = (-b - sqrt(b^2 - 8 * thetas[me]))/ 4. / thetas[me] 

	assert( p>= 0)
	return min(p, pmax)
end

#solve for equilibrium prices.  True model is log with a random, firm dependent effect
function solveNashPrices(thetas1::Array{Float64, 1}, thetas2::Array{Float64, 1}, 
						errors::Array{Float64, 2}, pmax; MAX_ITERS = 20, TOL=1e-6, trace=false, 
		prices=fill(convert(Float64, pmax / 2), 2))
	for i = 1:MAX_ITERS
		prices_old = copy(prices)
		prices[1] = bestResponsePrice(thetas1, prices_old, errors[1], pmax, :1)
		prices[2] = bestResponsePrice(thetas2, prices_old, errors[2], pmax, :2)
		dist = norm(prices - prices_old)
		if trace
			println(prices', dist)
		end
		if dist < TOL
			return prices
		end
	end
	error("Max Iters reached")
end

function margRev(thetas, prices, error, ID)
	me, you = indxID(ID)
	1./prices[me] + 2 * thetas[me] * prices[me] + thetas[you] * prices[you] + thetas[3] * error + thetas[4]
end

