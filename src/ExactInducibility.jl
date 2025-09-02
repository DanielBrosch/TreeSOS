using FlagSOS, JLD2, JuMP, LinearAlgebra, ProgressMeter, CairoMakie
using ClusteredLowRankSolver


function exact_certificate(T::BinaryTreeFlag, lvl::Int)
    T = labelCanonically(T)
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

t = BinaryTree
trees = generateAll(BinaryTreeFlag, 6)

cat4 = BinaryTreeFlag(t(t(false), t(t(false), t(t(false), t(false)))))
lvl = 4
(success, bound, sol, m) = exact_certificate(cat4, lvl)
@assert success && bound == 1
write_SOS_certificate(m, sol, cat4, lvl, "Inducibility_exact_Cat4.txt")

E4 = BinaryTreeFlag(t(t(t(false), t(false)), t(t(false), t(false))))
lvl = 4
(success, bound, sol, m) = exact_certificate(E4, lvl)
@assert success && bound == 3 // 7
write_SOS_certificate(m, sol, E4, lvl, "Inducibility_exact_E4.txt")

E5 = BinaryTreeFlag(t(t(t(false), t(false)), t(t(false), t(t(false), t(false)))))
lvl = 5
(success, bound, sol, m) = exact_certificate(E5, lvl)
@assert success && bound == 2 // 3
write_SOS_certificate(m, sol, E5, lvl, "Inducibility_exact_E5.txt")

E6 = BinaryTreeFlag(t(t(t(false), t(t(false), t(false))), t(t(false), t(t(false), t(false)))))
(success, bound, sol, m) = exact_certificate(E6, 6)
@assert success && bound == 10 // 31
write_SOS_certificate(m, sol, E6, 6, "Inducibility_exact_E6.txt")

if maxLvl >= 7
    E7 = BinaryTreeFlag(t(t(t(t(false), t(false)), t(t(false), t(false))), t(t(false), t(t(false), t(false)))))
    (success, bound, sol, m) = exact_certificate(E7, 7)
    @assert success && bound == 5 // 21
    write_SOS_certificate(m, sol, E7, 7, "Inducibility_exact_E7.txt")
end

if maxLvl >= 8
    E8 = BinaryTreeFlag(t(t(t(t(false), t(false)), t(t(false), t(false))), t(t(t(false), t(false)), t(t(false), t(false)))))
    (success, bound, sol, m) = exact_certificate(E8, 8)
    @assert success && bound == 45 // 889
    write_SOS_certificate(m, sol, E8, 8, "Inducibility_exact_E8.txt")
end

if maxLvl >= 9
    E9 = BinaryTreeFlag(labelCanonically(t(t(t(t(false), t(false)), t(t(false), t(false))), t(t(t(false), t(false)), t(t(false), t(t(false), t(false)))))))
    (success, bound, sol, m) = exact_certificate(E9, 9)
    @assert success && bound == 12 // 85
    write_SOS_certificate(m, sol, E9, 9, "Inducibility_exact_E9.txt")
end

if maxLvl >= 10
    E10 = BinaryTreeFlag(t(E5.tree, E5.tree))
    (success, bound, sol, m) = exact_certificate(E10, 10)
    @assert success && bound == 8 // 73
    write_SOS_certificate(m, sol, E10, 10, "Inducibility_exact_E10.txt")
end

## New exact bounds

if maxLvl >= 7
    T = BinaryTreeFlag(t(t(t(false), t(false)), cat4.tree))
    lvl = 7
    (success, bound, sol, m) = exact_certificate(T, lvl)
    @assert success && bound == 15 // 32
    write_SOS_certificate(m, sol, T, lvl, "Inducibility_exact_E2_Cat4.txt")
end

T = trees[15]
T == BinaryTreeFlag(t(t(t(false), t(false)), E4.tree))
lvl = 6
(success, bound, sol, m) = exact_certificate(T, lvl)
@assert success && bound == 45 // 217
write_SOS_certificate(m, sol, T, lvl, "Inducibility_exact_E2_E4.txt")