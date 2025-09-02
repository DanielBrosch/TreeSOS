## Experiments: Finding exact inducibility certificates and rounding them

using FlagSOS, JLD2, JuMP, LinearAlgebra, ProgressMeter, CairoMakie
using ClusteredLowRankSolver

## Load rounded bounds

# @load "Bounds/InducibilityBoundsRounded5.jld2" roundedBounds

bounds = Dict()

for lvl in 5:11
    @load "Bounds/InducibilityBoundsRounded$lvl.jld2" roundedBounds

    for (T, ineq) in roundedBounds
        T = BinaryTreeFlag(labelCanonically(T))
        if !haskey(bounds, T)
            bounds[T] = Vector{Any}()#Dict()
        end
        c = ineq.coeff[BinaryTreeFlag()]
        # bounds[T][lvl] = (Float64(c), c)
        push!(bounds[T], lvl => (Float64(c), c))
    end
end
boundsNotRounded = Dict()

for lvl in 5:11
    @load "Bounds/InducibilityBounds$lvl.jld2" inducibilityBounds

    for (T, c) in inducibilityBounds
        T = BinaryTreeFlag(labelCanonically(T))
        if !haskey(boundsNotRounded, T)
            boundsNotRounded[T] = Vector{Any}()#Dict()
        end
        # c = ineq.coeff[BinaryTreeFlag()]
        # bounds[T][lvl] = (Float64(c), c)
        push!(boundsNotRounded[T], lvl => (Float64(c), c))
    end
end



##

function exact_certificate(T::BinaryTreeFlag, lvl::Int)
    m = FlagModel{BinaryTreeFlag,:limit,Rational{Int}}()
    rM = addRazborovBlock!(m, lvl)
    computeSDP!(m)


    o = permute(BinaryTreeFlag(), 1:1)
    m.objective = -1 // 1 * T
    for _ in size(T):lvl-1
        m.objective *= 1 // 1 * o
    end

    obj, vars, blocks, blockSizes = buildStandardModel(m)

    varDicts = Dict(
        G => Dict{Any,Any}(
            mu[2] => Matrix(B[G])
            for (mu, B) in blocks if haskey(B, G) && blockSizes[mu] > 0
        )
        for G in vars
    )

    freeVarDicts = Dict(
        G => Dict(
            mu[2] => B[G][1, 1]
            for (mu, B) in blocks if haskey(B, G) && blockSizes[mu] < 0
        )
        for G in vars
    )

    @show freeVarDicts

    clObj = Objective(0, Dict(), Dict(:t => -1))
    clCons = Constraint[]

    for G in vars
        c = Constraint(get(obj.coeff, G, 0), varDicts[G], Dict(:t => 1))
        push!(clCons, c)
    end
    P = Problem(Minimize(clObj), clCons)

    bits = 2^10
    setprecision(bits)
    status, primalsol, dualsol, time, errorcode = solvesdp(P;
        # step_length_threshold=BigFloat(10e-14),
        duality_gap_threshold=1 / (BigFloat(2)^(bits / 4)))
    objvalue(P, dualsol)

    success, exactdualsol = exact_solution(P, primalsol, dualsol; settings=RoundingSettings(kernel_lll=false, kernel_use_primal=false))
    # success, exactdualsol, transform = exact_solution(P, primalsol, dualsol; settings = RoundingSettings(kernel_lll = false, kernel_use_primal=false), transformed=true)

    @assert success
    success, objvalue(P, exactdualsol), exactdualsol, m
end

function write_SOS_certificate(m, sol, T, lvl, file=nothing)
    res = Dict{Any,Any}(T => Rational{BigInt}.(M) for (T, M) in sol.matrixvars)
    ct = Rational{BigInt}(sol.freevars[:t])
    quot = m.subModels[1].quotient
    for i in eachindex(quot)
        res["Q$i"] = ct
    end
    if lvl > size(T)
        o = permute(BinaryTreeFlag(), 1:1)
        Ti = -1 // 1 * T
        for _ in size(T):lvl-1
            for (t, c) in Ti.coeff
                to = -1 // 1 * t + t * o
                k = findfirst(x -> x == to, quot)
                res["Q$k"] -= c
            end
            Ti *= 1 // 1 * o
        end
    end
    res
    m.objective = -1 // 1 * T
    if file === nothing
        return FlagSOS.verifySOS(m, [res])
    else
        !isdir("Certificates") && mkdir("Certificates")
        filepath = joinpath("Certificates", file)
        open(filepath, "w") do io
            return FlagSOS.verifySOS(m, [res]; io=io)
        end
    end
end

##

function even_tree_inducibility(n)
    function c(n)
        n == 0 && return 1
        n == 1 && return 1
        if n % 2 == 0
            s = Int(n / 2)
            return (c(s)^2) // (2^(2 * s) - 2)
        else
            s = Int((n - 1) / 2)
            return (c(s) * c(s + 1)) // (2^(2 * s) - 1)
        end
    end
    return c(n) * factorial(n)
end

##

trees = generateAll(BinaryTreeFlag, 6)

cat4 = trees[5]
bounds[cat4]
lvl = 4
(success, bound, sol, m) = exact_certificate(cat4, lvl)
write_SOS_certificate(m, sol, cat4, lvl, "Inducibility_Cat4.txt")

E4 = trees[6]
lvl = 4
bounds[E4]
(success, bound, sol, m) = exact_certificate(E4, lvl)
write_SOS_certificate(m, sol, E4, lvl)

E5 = trees[8]
bounds[E5]
boundsNotRounded[E5]
lvl = 5
(success, bound, sol, m) = exact_certificate(E5, lvl)
write_SOS_certificate(m, sol, E5, lvl)



E6 = trees[14]
bounds[E6]
exact_certificate(E6, 6)



t = BinaryTree
E7 = BinaryTreeFlag(labelCanonically(t(t(t(false), t(t(false), t(false))), t(t(t(false), t(false)), t(t(false), t(false))))))
bounds[E7]
exact_certificate(E7, 7)

E8 = BinaryTreeFlag(labelCanonically(t(t(t(t(false), t(false)), t(t(false), t(false))), t(t(t(false), t(false)), t(t(false), t(false))))))
bounds[E8]
exact_certificate(E8, 8)

E9 = BinaryTreeFlag(labelCanonically(t(t(t(t(false), t(false)), t(t(false), t(false))), t(t(t(false), t(false)), t(t(false), t(t(false), t(false)))))))
bounds[E9]
exact_certificate(E9, 9)

E10 = BinaryTreeFlag(labelCanonically(t(E5.tree, E5.tree)))
bounds[E10]
exact_certificate(E10, 10)

## New exact bounds

T = trees[11]
bounds[T]
lvl = 7
(success, bound, sol, m) = exact_certificate(T, lvl)
write_SOS_certificate(m, sol, T, lvl)

T = trees[15]
bounds[T]
lvl = 6
(success, bound, sol, m) = exact_certificate(T, lvl)
write_SOS_certificate(m, sol, T, lvl)

## Problem cases

T = trees[12]
bounds[T]
exact_certificate(T, 8) # not exact? fails at lvl 8 still

T = trees[13]
bounds[T]
exact_certificate(T, 8) # can't round at lvl 8, not exact?, fails at lvl 8.

TIrr = trees[9]
bounds[TIrr]
lvl = 10
(success, bound, sol, m) = exact_certificate(TIrr, lvl)
write_SOS_certificate(m, sol, TIrr, lvl)
##

# E5 = trees[8] has a phantom edge?
# irrational tree = trees[9]
trees = generateAll(BinaryTreeFlag, 5)[8:8]#[5:5]#[end:end]
trees = labelCanonically.(trees)

## Fractalizer checks
using Combinatorics

function replace_leaves(T::BinaryTree, replacements::Vector{BinaryTree})
    @assert size(T) == length(replacements)

    if size(T) == 1
        # T.left = replacements[1].left
        # T.right = replacements[1].right
        return replacements[1]
    else
        k = size(T.left)
        TLeft = replace_leaves(T.left, replacements[1:k])
        TRight = replace_leaves(T.right, replacements[k+1:end])
        return BinaryTree(TLeft, TRight)
    end
end

function subtrees_in_fractal(T::BinaryTree, n::Int)
    res = Dict{BinaryTree,Rational{Int}}()
    if n == 0
        return Dict(BinaryTree(true) => 1 // 1)
    elseif n == 1
        return Dict(BinaryTree(false) => 1 // 1)
    elseif n == 2
        return Dict(BinaryTree(BinaryTree(false), BinaryTree(false)) => 1 // 1)
    end
    multicomb = with_replacement_combinations(1:size(T), n)
    k = length(multicomb)
    # fact = k * (1 - 1 // (size(T)^(n-1)))
    # fact = k * (1 - size(T) // k)
    # fact = k * (1 - size(T) // (size(T)^(n)))
    fact = size(T)^n * (1 - 1 // (size(T)^(n - 1)))
    for c in multicomb
        # @show c
        counts = [count(c .== i) for i in 1:size(T)]
        # @show counts
        T_sub = FlagSOS.subFlag(T, [i for i in 1:size(T) if counts[i] != 0])
        # @show T_sub
        filter!(!iszero, counts)
        # @show counts
        length(counts) == 1 && continue

        fact2 = factorial(n) // prod(factorial.(counts))

        sub_counts = [subtrees_in_fractal(T, c) for c in counts]
        for leaf_replacements in Iterators.product(sub_counts...)
            c = prod(ci for (_, ci) in leaf_replacements)
            trees = [t for (t, _) in leaf_replacements]
            new_tree = labelCanonically(replace_leaves(T_sub, trees))
            res[new_tree] = get(res, new_tree, 0 // 1) + fact2 * c // fact
        end
    end
    # @show T, n
    # @show res
    # @show fact
    # @show k
    @assert sum(values(res)) == 1
    # @show collect(multicomb)
    return res
end



res = subtrees_in_fractal(trees[3].tree, 6)
res[trees[15].tree] # sharp!

res = subtrees_in_double_cat(1//2, 6)
res[trees[11].tree] # sharp

res = subtrees_in_fractal(trees[9].tree, 5)
res[trees[9].tree]
##

function caterpillar(n::Int)
    if n == 0
        return BinaryTree(true)
    elseif n == 1
        return BinaryTree(false)
    else
        return BinaryTree(BinaryTree(false), caterpillar(n - 1))
    end
end

# Infinite tree: Left of root: E_Infty, right of root: Cat_infty. 
function subtrees_in_EInfty_CatInfty(E_ratio::Rational{Int}, n::Int)
    res = Dict{BinaryTree,Rational{Int}}()
    E2 = BinaryTree(BinaryTree(false), BinaryTree(false))
    for no_left_verts in 0:n
        prob = binomial(n, no_left_verts) * E_ratio^no_left_verts * (1 - E_ratio)^(n - no_left_verts)

        left_trees = subtrees_in_fractal(E2, no_left_verts)
        for (T, c) in left_trees
            if no_left_verts == n
                new_tree = T
            elseif no_left_verts == 0
                new_tree = caterpillar(n)
            else
                new_tree = labelCanonically(BinaryTree(T, caterpillar(n - no_left_verts)))
            end
            @show no_left_verts, T, new_tree
            res[new_tree] = get(res, new_tree, 0 // 1) + prob * c
        end
    end
    return res
end

function subtrees_in_EInfty_EInfty(E_ratio::Rational{Int}, n::Int)
    res = Dict{BinaryTree,Rational{Int}}()
    E2 = BinaryTree(BinaryTree(false), BinaryTree(false))
    for no_left_verts in 0:n
        prob = binomial(n, no_left_verts) * E_ratio^no_left_verts * (1 - E_ratio)^(n - no_left_verts)

        left_trees = subtrees_in_fractal(E2, no_left_verts)
        right_trees = subtrees_in_fractal(E2, n - no_left_verts)
        for ((S, d),(T, c)) in Iterators.product(left_trees, right_trees)
            if no_left_verts == n
                new_tree = S
            elseif no_left_verts == 0
                new_tree = T
            else
                new_tree = labelCanonically(BinaryTree(S,T))
            end
            @show no_left_verts, T, new_tree
            res[new_tree] = get(res, new_tree, 0 // 1) + prob * c*d
        end
    end
    return res
end

function subtrees_in_double_cat(ratio::Rational{Int}, n::Int)
    res = Dict{BinaryTree,Rational{Int}}()
    for no_left_verts in 0:n
        prob = binomial(n, no_left_verts) * ratio^no_left_verts * (1 - ratio)^(n - no_left_verts)

        if no_left_verts == n
            new_tree = caterpillar(n)
        elseif no_left_verts == 0
            new_tree = caterpillar(n)
        else
            new_tree = labelCanonically(BinaryTree(caterpillar(no_left_verts), caterpillar(n - no_left_verts)))
        end
        res[new_tree] = get(res, new_tree, 0 // 1) + prob
    end
    return res
end


tmp = [subtrees_in_EInfty_CatInfty(p, 6)[trees[11].tree] for p in 1//100:1//100:99//100]
tmp = [subtrees_in_EInfty_CatInfty(p, 5)[trees[9].tree] for p in 1//100:1//100:99//100]
tmp = [subtrees_in_EInfty_EInfty(p, 5)[trees[9].tree] for p in 1//100:1//100:99//100]
tmp = [subtrees_in_double_cat(p, 6)[trees[11].tree] for p in 1//100:1//100:99//100]
# tmp = [subtrees_in_double_cat(p, 5)[trees[9].tree] for p in 1//100:1//100:99//100]

##




##

mus = collect(keys(blocks))

mu = mus[3]

FF = find_field(primalsol, dualsol, 2)

# full kernel
full_kernel = ClusteredLowRankSolver.detecteigenvectors(primalsol.matrixvars[mu], dualsol.matrixvars[mu]; settings=RoundingSettings(kernel_use_primal=false))
# expected kernel
expected_kernel = ClusteredLowRankSolver.detecteigenvectors(primalsol.matrixvars[mu], dualsol.matrixvars[mu]; settings=RoundingSettings(kernel_use_primal=true))

expected_kernel[1] += 2 * expected_kernel[2]
full_kernel[1] += full_kernel[2]
full_kernel[4] += full_kernel[2]
full_kernel[1] += full_kernel[4]

expected_kernel
full_kernel


for b in expected_kernel
    basisEl = sum(c * G for (c, G) in zip(Rational{Int}.(b), rM.basis[mu[2]]))
    display(basisEl)
end

for b in full_kernel
    basisEl = sum(c * G for (c, G) in zip(Rational{Int}.(b), rM.basis[mu[2]]))
    display(basisEl)
end

##

success, exactdualsol = exact_solution(P, primalsol, dualsol; settings=RoundingSettings(kernel_lll=true, kernel_errbound=1e-8, kernel_round_errbound=1e-6, kernel_bits=256, kernel_use_primal=false), FF=FF[1])
success, exactdualsol = exact_solution(P, primalsol, dualsol; settings=RoundingSettings(kernel_lll=false, kernel_errbound=1e-8, kernel_round_errbound=1e-6, kernel_bits=256, kernel_use_primal=false), FF=FF[1])
success, exactdualsol = exact_solution(P, primalsol, dualsol; settings=RoundingSettings(kernel_lll=true, kernel_errbound=1e-8, kernel_round_errbound=1e-6, kernel_bits=256, kernel_use_primal=false))
success, exactdualsol = exact_solution(P, primalsol, dualsol; settings=RoundingSettings(kernel_lll=true, kernel_errbound=1e-8, kernel_round_errbound=1e-6, kernel_bits=256, kernel_use_primal=false), transformed=true)
success, exactdualsol = exact_solution(P, primalsol, dualsol; settings=RoundingSettings(kernel_lll=true, kernel_errbound=1e-8, kernel_round_errbound=1e-6, kernel_bits=256, kernel_use_primal=false))
success, exactdualsol = exact_solution(P, primalsol, dualsol; settings=RoundingSettings(kernel_lll=false, kernel_errbound=1e-6, kernel_round_errbound=1e-8, kernel_bits=128, kernel_use_primal=true))
success, exactdualsol = exact_solution(P, primalsol, dualsol; settings=RoundingSettings(kernel_lll=true, kernel_use_primal=false))
success, exactdualsol = exact_solution(P, primalsol, dualsol; settings=RoundingSettings(kernel_lll=true, kernel_use_primal=true))
success, exactdualsol = exact_solution(P, primalsol, dualsol; settings=RoundingSettings(kernel_lll=false, kernel_use_primal=true))
success, exactdualsol = exact_solution(P, primalsol, dualsol; settings=RoundingSettings(kernel_lll=false, kernel_use_primal=true), transformed=true)
success, exactdualsol = exact_solution(P, primalsol, dualsol; settings=RoundingSettings(kernel_lll=false, kernel_use_primal=true))
success, exactdualsol = exact_solution(P, primalsol, dualsol; settings=RoundingSettings(kernel_lll=true, kernel_use_primal=true), FF=FF[1])

success, exactdualsol = exact_solution(P, primalsol, dualsol; settings=RoundingSettings(kernel_lll=false, kernel_use_primal=false))
# success, exactdualsol, transform = exact_solution(P, primalsol, dualsol; settings = RoundingSettings(kernel_lll = false, kernel_use_primal=false), transformed=true)

@assert success
objvalue(P, exactdualsol)

##

for (mu, M) in exactdualsol.matrixvars
    println()
    @show mu
    display(M)
    # display(primalsol.matrixvars[mu])
    # display(eigen(primalsol.matrixvars[mu]))
    display(dualsol.matrixvars[mu])
    # display(eigen(BigFloat.(M)))
    # println("In basis:")
    # display(transform[mu])
    # for b in 1:size(transform[mu],2)
    #     basisEl = sum(c*G for (c,G) in zip(Rational{Int}.(transform[mu][:, b]),rM.basis[mu[2]]))
    #     display(basisEl)
    # end
end

for (mu, M) in exactdualsol.freevars
    println()
    @show mu
    display(M)
end

##

x1vec = [1, -2]
mu = BinaryTreeFlag(false)

for (T, B) in rM.sdpData
    c = dot(B[mu], x1vec * x1vec')
    @show T, c
end

x2vec = [1, -1, 0, 1, 1]
mu = BinaryTreeFlag(BinaryTree(BinaryTree(false), BinaryTree(BinaryTree(false), BinaryTree(false))))

for (T, B) in rM.sdpData
    c = dot(B[mu], x2vec * x2vec')
    @show T, c
end

##
using Hypatia
m.objective = -1 * T
jm = buildJuMPModel(m)

set_optimizer(jm.model, Hypatia.Optimizer)
optimize!(jm.model)
@show T
# @show termination_status(jm.model)
@show objective_value(jm.model)

# using FlagSOS, JLD2, JuMP, LinearAlgebra, ProgressMeter, CairoMakie
# using ClusteredLowRankSolver

# # E5 = trees[8] has a phantom edge?
# # irrational tree = trees[9]
# trees = generateAll(BinaryTreeFlag, 5)[8:8]#[5:5]#[end:end]
# trees = labelCanonically.(trees)


# m = FlagModel{BinaryTreeFlag,:limit,Rational{Int}}()

# lvl = 5

# rM = addRazborovBlock!(m, lvl)
# computeSDP!(m)


# modelSize(m) # lvl 9: 70₁54₃13₁₁9₁1₄₆

# #


# T = trees[1]
# o = permute(BinaryTreeFlag(),1:1)
# m.objective = -1//1 * T
# for _ in size(T):lvl-1 
#     m.objective *= 1//1*o 
# end


# obj, vars, blocks, blockSizes = buildStandardModel(m)
# # constant_one = 1//1*o 
# # for _ = 1:lvl-1 
# #     constant_one *= 1//1*o 
# # end

# # constant_one

# varDicts = Dict(
#     G => Dict{Any,Any}(
#         mu => Matrix(B[G])
#         for (mu, B) in blocks if haskey(B, G) && blockSizes[mu] > 0
#     )
#     for G in vars
# )

# # bound variables
# # for G in vars
# #     if G != o
# #         varDicts[G][(:lower, G)] = [1;;]
# #         varDicts[G][(:upper, G)] = [-1;;]
# #         varDicts[o][(:upper, G)] = [1;;]
# #     end
# # end

# freeVarDicts = Dict(
#     G => Dict(
#         mu => B[G][1, 1]
#         for (mu, B) in blocks if haskey(B, G) && blockSizes[mu] < 0
#     )
#     for G in vars
# )



# clObj = Objective(0, Dict(), Dict(:t => -1))
# clCons = Constraint[]

# for G in vars
#     c = Constraint(get(obj.coeff, G, 0), varDicts[G], Dict(:t=>1))
#     # @show c
#     push!(clCons, c)
# end
# P = Problem(Minimize(clObj), clCons)

# bits = 2^9
# setprecision(bits)
# status, primalsol, dualsol, time, errorcode = solvesdp(P;
#     # step_length_threshold=BigFloat(10e-14),
#     duality_gap_threshold=1 / (BigFloat(2)^(bits / 4)));
# objvalue(P, dualsol)

# ##

# mus = collect(keys(blocks))

# mu = mus[3]

# FF = find_field(primalsol, dualsol, 2)

# # full kernel
# full_kernel = ClusteredLowRankSolver.detecteigenvectors(primalsol.matrixvars[mu], dualsol.matrixvars[mu]; settings = RoundingSettings(kernel_use_primal=false))
# # expected kernel
# expected_kernel = ClusteredLowRankSolver.detecteigenvectors(primalsol.matrixvars[mu], dualsol.matrixvars[mu]; settings = RoundingSettings(kernel_use_primal=true))

# expected_kernel[1] += 2*expected_kernel[2]
# full_kernel[1] += full_kernel[2]
# full_kernel[4] += full_kernel[2]
# full_kernel[1] += full_kernel[4]

# expected_kernel
# full_kernel


# for b in expected_kernel
#     basisEl = sum(c*G for (c,G) in zip(Rational{Int}.(b),rM.basis[mu[2]]))
#     display(basisEl)
# end

# for b in full_kernel
#     basisEl = sum(c*G for (c,G) in zip(Rational{Int}.(b),rM.basis[mu[2]]))
#     display(basisEl)
# end

# ##

# success, exactdualsol = exact_solution(P, primalsol, dualsol; settings = RoundingSettings(kernel_lll = true, kernel_errbound=1e-8,kernel_round_errbound=1e-6, kernel_bits = 256, kernel_use_primal=false), FF = FF[1])
# success, exactdualsol = exact_solution(P, primalsol, dualsol; settings = RoundingSettings(kernel_lll = false, kernel_errbound=1e-8,kernel_round_errbound=1e-6, kernel_bits = 256, kernel_use_primal=false), FF = FF[1])
# success, exactdualsol = exact_solution(P, primalsol, dualsol; settings = RoundingSettings(kernel_lll = true, kernel_errbound=1e-8,kernel_round_errbound=1e-6, kernel_bits = 256, kernel_use_primal=false))
# success, exactdualsol = exact_solution(P, primalsol, dualsol; settings = RoundingSettings(kernel_lll = true, kernel_errbound=1e-8,kernel_round_errbound=1e-6, kernel_bits = 256, kernel_use_primal=false), transformed=true)
# success, exactdualsol = exact_solution(P, primalsol, dualsol; settings = RoundingSettings(kernel_lll = true, kernel_errbound=1e-8,kernel_round_errbound=1e-6, kernel_bits = 256, kernel_use_primal=false))
# success, exactdualsol = exact_solution(P, primalsol, dualsol; settings = RoundingSettings(kernel_lll = false, kernel_errbound=1e-6,kernel_round_errbound=1e-8, kernel_bits = 128, kernel_use_primal=true))
# success, exactdualsol = exact_solution(P, primalsol, dualsol; settings = RoundingSettings(kernel_lll = true, kernel_use_primal=false))
# success, exactdualsol = exact_solution(P, primalsol, dualsol; settings = RoundingSettings(kernel_lll = true, kernel_use_primal=true))
# success, exactdualsol = exact_solution(P, primalsol, dualsol; settings = RoundingSettings(kernel_lll = false, kernel_use_primal=true))
# success, exactdualsol = exact_solution(P, primalsol, dualsol; settings = RoundingSettings(kernel_lll = false, kernel_use_primal=true), transformed=true)
# success, exactdualsol = exact_solution(P, primalsol, dualsol; settings = RoundingSettings(kernel_lll = false, kernel_use_primal=true))
# success, exactdualsol = exact_solution(P, primalsol, dualsol; settings = RoundingSettings(kernel_lll = true, kernel_use_primal=true), FF = FF[1])

# success, exactdualsol = exact_solution(P, primalsol, dualsol; settings = RoundingSettings(kernel_lll = false, kernel_use_primal=false))
# # success, exactdualsol, transform = exact_solution(P, primalsol, dualsol; settings = RoundingSettings(kernel_lll = false, kernel_use_primal=false), transformed=true)

# @assert success
# objvalue(P, exactdualsol)

# ##

# for (mu, M) in exactdualsol.matrixvars
#     println()
#     @show mu 
#     display(M)
#     # display(primalsol.matrixvars[mu])
#     # display(eigen(primalsol.matrixvars[mu]))
#     display(dualsol.matrixvars[mu])
#     # display(eigen(BigFloat.(M)))
#     # println("In basis:")
#     # display(transform[mu])
#     # for b in 1:size(transform[mu],2)
#     #     basisEl = sum(c*G for (c,G) in zip(Rational{Int}.(transform[mu][:, b]),rM.basis[mu[2]]))
#     #     display(basisEl)
#     # end
# end

# for (mu, M) in exactdualsol.freevars 
#     println()
#     @show mu 
#     display(M)
# end

# ##

# x1vec = [1, -2]
# mu = BinaryTreeFlag(false)

# for (T, B) in rM.sdpData 
#     c = dot(B[mu], x1vec*x1vec')
#     @show T, c
# end

# x2vec = [1, -1, 0, 1, 1]
# mu = BinaryTreeFlag(BinaryTree(BinaryTree(false), BinaryTree(BinaryTree(false), BinaryTree(false))))

# for (T, B) in rM.sdpData 
#     c = dot(B[mu], x2vec*x2vec')
#     @show T, c
# end

# ##
# using Hypatia
# m.objective = -1 * T
# jm = buildJuMPModel(m)

# set_optimizer(jm.model, Hypatia.Optimizer)
# optimize!(jm.model)
# @show T
# # @show termination_status(jm.model)
# @show objective_value(jm.model)