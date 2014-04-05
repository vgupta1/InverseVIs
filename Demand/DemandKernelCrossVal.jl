using DataFrames

#cd("/Users/VGupta/Documents/Research/Julia Stuff/Demand")
include("BertrandNash2.jl");  # does the generation of the prices
include("fitNash2.jl");  # does the fitting

thetas1 = [-1.2, .5, 1, -10]
thetas2 = [.3, -1, 1, -10]
pmax = .45; N = 200;

function genErrors( ; mu1 = 5, mu2 = 5, mu3 = 5, sigma1 = 1.5, sigma2 = 1.5, sigma3 = 1.5)
    eps1 = rand(Normal(mu1, sigma1)); eps2 = rand(Normal(mu2, sigma2)); eps3 = rand(Normal(mu3, sigma3))
    eps3 = (eps1 + eps2 + eps3)/3
	return eps1, eps2, eps3
end

srand(8675309)
inPrices = Array(Float64, N, 2)
inErrors = Array(Float64, N, 3)
for iRun = 1:N
    inErrors[iRun, :] = [genErrors()...]
    pstar = solveNashPrices(thetas1, thetas2, inErrors[iRun, :], pmax, prices=[.3, .3])
	inPrices[iRun, :] = pstar'
end
indx = find( abs(inErrors[:, 3] - median(inErrors[:, 3]) ) .<= 1e-2 )[1];

function fit(inPrices, inErrors, c, lambda_reg)
    ### Thin out the grids
    # p1_grid = sort(inPrices[:, 1]); p2_grid = sort(inPrices[:, 2])
    # p1_grid = thin(p1_grid, .45/50, pmax)
    # p2_grid = thin(p2_grid, .45/50, pmax)
    p1_grid = p2_grid = linspace(.1, pmax, 50)
    
    #Augment the data as necessary
    N1 = length(p1_grid); N2 = length(p2_grid); totN = N + N1 + N2
    data = [inPrices inErrors[:, 3]]
    data = [data; p1_grid repmat([inPrices[indx, 2]], N1, 1) repmat([inErrors[indx, 3]], N1, 1)]
    data = [data; repmat([inPrices[indx, 1]], N2, 1) p2_grid repmat([inErrors[indx, 3]], N2, 1)]
    data = [data; p1_grid[1] inPrices[indx, 2] inErrors[indx, 3];
                  inPrices[indx, 1] p2_grid[1] inErrors[indx, 3]]
 
    #build a kernel, the graham matrix, and it's cholesky
    #only use significant eigenvalues
    k(x,y) = exp( -1/c * norm(x-y)^2 )
    K = [k(data[i, :], data[j, :]) for i =1:N+N1+N2+2, j=1:totN+2]
    D, V = eig(K)
    tailN = find( 1 - (cumsum(D) / sum(D)) .<= 1-1e-6)[1]
    Dbar = copy(D)
    Dbar[1:tailN-1] = 0;
    Vbar = V * diagm(sqrt(Dbar));
    
    m, vZ1, vZ2, vf1, vf2, vIntercept = setUpFitting(totN+2, Vbar[:, tailN:totN+2])

    normalize(m, vf1[totN + 1], margRev(thetas1, data[totN + 1, 1:2], data[totN + 1, 3], :1))                
    normalize(m, vf2[totN + 2], margRev(thetas2, data[totN + 2, 1:2], data[totN + 2, 3], :2))                
    addDecreasing(m, vf1, (N+1):N+N1, TOL=1e-6)
    addDecreasing(m, vf2, (N+N1+1):totN, TOL=1e-6)
    
    resids = Variable[]
    @defVar(m, maxVal >=0)
    for iRun = 1:N
        resid = addResids(m, vf1, vf2, iRun, data, pmax)
        push!(resids, resid)
        @addConstraint(m, maxVal >= resid)
    end

    @defVar(m, reg_term[1:2])
    HNorm1 = QuadExpr(vZ1[:], vZ1[:], ones(length(vZ1)), AffExpr())
    HNorm2 = QuadExpr(vZ2[:], vZ2[:], ones(length(vZ2)), AffExpr())
    addConstraint(m, HNorm1 <= reg_term[1]) 
    addConstraint(m, HNorm2 <= reg_term[2])

    @setObjective(m, Min, sum{resids[i], i=1:N}/N + lambda_reg * (reg_term[1] + reg_term[2]))
    println(solve(m))

    #extract solution
    fZ1 = [getValue(vZ1[i])::Float64 for i =1:length(vZ1)]
    fZ2 = [getValue(vZ2[i])::Float64 for i = 1:length(vZ2)]
    fIntercept = [getValue(vIntercept[i])::Float64 for i=1:2]
    fMR1_fun = createfun(fZ1, fIntercept[1], Vbar[:, tailN:totN+2], data, k)
    fMR2_fun = createfun(fZ2, fIntercept[2], Vbar[:, tailN:totN+2], data, k)

    return fMR1_fun, fMR2_fun
end

srand(516746226)
function testOut(fMR1_fun, fMR2_fun, Nout, METHOD)
    if METHOD == :pdiffs
        fPrices = Array(Float64, Nout, 2)
        oPrices = Array(Float64, Nout, 2)
        oErrors = Array(Float64, Nout, 3)
        for iRun = 1:Nout
            oErrors[iRun, :] = [genErrors()...]
            pstar = solveNashPrices(thetas1, thetas2, oErrors[iRun, :], pmax)
            oPrices[iRun, :] = pstar'
        
           println(iRun)
            
            fPrices[iRun, :] = solveNashPrices(fMR1_fun, fMR2_fun, oErrors[iRun, 3], pmax, 
                                                   prices=copy(pstar), trace=false)
        end
        price_diffs = oPrices - fPrices
        return norm(price_diffs) / sqrt(Nout)
    else
        eps = 0.
        oErrors = Array(Float64, Nout, 3)
        for iRun = 1:Nout
            oErrors[iRun, :] = [genErrors()...]
            pstar = solveNashPrices(thetas1, thetas2, oErrors[iRun, :], pmax)
            MR1 = fMR1_fun(pstar..., oErrors[iRun, 3])
            MR2 = fMR2_fun(pstar..., oErrors[iRun, 3])

            y1 = max(0, MR1); y2 = max(0, MR2)
            eps += max(0, pmax * (y1 + y2) - pstar[1] * MR1 -pstar[2] * MR2 )
        end
        return eps/Nout
    end
end

lamb_grid = 10. .^(-6:.5:1)
c_grid = 10. .^(-3:.25:1)

#1 REdo this with tighter lambda, larger c
# lamb_grid = 10. .^(-6:.5:-2)
# c_grid = 10. .^(1:.5:3)
fold_dict = Dict{(Float64, Float64), Array{Float64, 1}}()

NSplits = 10
for (ix, l) in enumerate(lamb_grid)
    for (jx, c) in enumerate(c_grid)
        if haskey(fold_dict, (l,c))
            continue
        end
        
        fMR1_fun, fMR2_fun = fit(inPrices, inErrors, c, l)
        errors = zeros(NSplits)
        good_indices = Int[]
        for iSplit = 1:NSplits
            try
                errors[iSplit] = testOut(fMR1_fun, fMR2_fun, 200, :pdiffs)
                push!(good_indices, iSplit)
            catch y
                println("Error for $l $c \t", y)
            end
        end
        println("Good_indices: $(length(good_indices))")
        println(good_indices)
        println(errors)
        println(errors[good_indices])
        if length(good_indices) > 0 
            s = length(good_indices) > 1 ? std(errors[good_indices]) : -1.
            fold_dict[(l, c)] = [mean(errors[good_indices]), s]
        end
    end
end

#Create an array of the relevant values
fold_vals = Array(Float64, 0, 4)
for l = lamb_grid
    for c = c_grid
        if haskey(fold_dict, (l, c))
            fold_vals = [fold_vals ; [l c fold_dict[(l,c)][1] fold_dict[(l,c)][2] ] ]
            println("$l \t $c \t $(fold_dict[(l,c)][1]) \t $(fold_dict[(l,c)][2])")
        end
    end
end

df = DataFrame(fold_vals)
colnames!(df, ["Lambda", "C", "Mean", "SD"])
writetable("crossVal_pdiffs.csv", df)

