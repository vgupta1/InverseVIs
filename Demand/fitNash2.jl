## File contains the fitting funcs with Kernel


using JuMP
using Gurobi
using Distributions
using Roots

#size(Vbar) = N x 40
#Vbar Vbar' approx K
#Vbar' alpha = z
#Vbar z = f

function createfun(fZ, fintercept, Vbar, data, k )
	#construct the alphas
	alphas = Vbar' \ fZ
	println(size(alphas))
	(p1, p2, gdp) -> (sum([alphas[i]*k(data[i, :], [p1 p2 gdp]) for i=1:length(alphas)]) + fintercept)
end 
	

function setUpFitting(N, Vbar)
	m = Model(solver=GurobiSolver(OutputFlag=false))
	@defVar(m, vZ1[1:size(Vbar, 2)])
	@defVar(m, vZ2[1:size(Vbar, 2)])
	@defVar(m, vf1[1:N])
	@defVar(m, vf2[1:N])
	@defVar(m, vIntercept[1:2])

	for i = 1:N
		@addConstraint(m, vf1[i] == sum{Vbar[i, j]*vZ1[j], j=1:size(Vbar, 2)} + vIntercept[1])
		@addConstraint(m, vf2[i] == sum{Vbar[i, j]*vZ2[j], j=1:size(Vbar, 2)} + vIntercept[2])
	end
	return m, vZ1, vZ2, vf1, vf2, vIntercept
end

#insist that the MR at the mean point should be equal to gien  quantity for each
function normalize(m, vfj, val)
	@addConstraint(m, vfj  == val)
end

#indices correspond to the section of data where relevant price is changing, otehrs fixed to median
function addDecreasing(m, vf, indices; TOL =0)
	for i = indices[1:length(indices) - 1]
		@addConstraint(m, vf[i] >= vf[i + 1] - TOL)
	end
	#VG: Should add one more for pmax
end

function addPositivity(m, vf, indx; TOL=0)
	@addConstraint(m, vf[indx] >= TOL)
end


#adds the constraints to the model for the residual for this observation
function addResids(m, vf1, vf2, indx, price_data, pmax)
	@defVar(m, resid >= 0)
	@defVar(m, y[1:2] >= 0)
	@addConstraint(m, y[1] >= vf1[indx])
	@addConstraint(m, y[2] >= vf2[indx])
	@addConstraint(m, pmax * (y[1] + y[2]) 
					- price_data[indx, 1] * vf1[indx] - price_data[indx, 2] * vf2[indx] <= resid)
	return resid
end

#calc an out of sample residual.
function calcResid(fMR1, fMR2, data, pmax)
	y = Array(Float64, 2)
	y[1] = max(0, fMR1(data...))
	y[2] = max(0, fMR2(data...))
	return (pmax * (y[1] + y[2]) - data[1] * fMR1(data...) - data[2] * fMR2(data...))
end

#VG Figure out how to make it call this one instead....
function bestResponsePrice2(fMR, p1, p2, error, pmax, ID)
	if ID == :1
		pstar = customRoot(p-> fMR(p, p2, error), pmax, p1)  #works really well with exp features
	else	
		pstar = customRoot(p-> fMR(p1, p, error), pmax, p2)   #works really well with exp features
	end
	# assert( pstar >= 0 )
	pstar = max(0, pstar)
	return min(pstar, pmax)
end

#solve for equilibrium prices.  True model is log with a random, firm dependent effect
function solveNashPrices(fMR1, fMR2, error, pmax; MAX_ITERS = 20, TOL=1e-6, trace=false, 
		prices=fill(0.5 * pmax, 2) )
	dist = 0
	for i = 1:MAX_ITERS
		prices_old = copy(prices)
		prices[1] = bestResponsePrice2(fMR1, prices_old..., error, pmax, :1)
		prices[2] = bestResponsePrice2(fMR2, prices_old..., error, pmax, :2)
		dist = norm(prices - prices_old)
		if trace
			println([prices; dist]')
		end
		if dist < TOL
			return prices
		end
	end
	warn("Max Iters Reached \t $dist")
	return prices
end

# WTF.  How do you optimize gaussian functions?
# function customRoot(f, pmax; NUMPTS=500)
# 	#First attempt, will return "some" zero, multiplicity is an issue
# 	if f(1e-3) * f(pmax) < 0
# 		return fzero(f, [1e-3, pmax])
# 	end

# 	#sample a bunch of points on (0, pmax)
# 	ps = [1e-3; sort(pmax * rand(NUMPTS)); pmax]
# 	assert(ps[1] > 0)
# 	#prefer low-valued prices
# 	if minimum(ps) > 0
# 		warn("Potential Exhaustive Search Failure")	
# 		return pmax
# 	else
# 		indx = find( ps .< 0 )[1]
# 		assert(ps[indx-1] > 0)
# 		return fzero(f, ps[(indx-1):indx])
# 	end
# end


#assumes f(0) > 0 and we have decreasingness
function customRoot(f, pmax, startp; MAX_BRACKET = 100)
	if f(pmax) > 0
		return pmax
	end
	if f(startp) <= 0
		return fzero(f, [.1, startp])
	else
		return fzero(f, [startp, pmax]) 
	end
end

#Try thinning out p1_grid, p2_grid
function thin(p1_grid, delta, pmax)
    p1_gridp = Float64[]; push!(p1_gridp, p1_grid[1])
    for i = 2:N
        if p1_grid[i] - p1_gridp[length(p1_gridp)] > delta
            push!(p1_gridp, p1_grid[i] )
        end
    end
    if abs(p1_gridp[length(p1_gridp)] - pmax) > .05
        push!(p1_gridp, pmax)
    end
    return p1_gridp
end


