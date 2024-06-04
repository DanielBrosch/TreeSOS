using FlagSOS, JLD2, JuMP, LinearAlgebra, ProgressMeter, CairoMakie

## Computes outer approximations of tree-profiles

function boundProfileQuantumFlag(m::FlagModel{T,N,D}, G1::QuantumFlag, G2::QuantumFlag, razLevel; steps=range(0 // 1, 1 // 1, 10), solver=solver, type=Float64, lb=true, ub=true, coeffLimit=-1, maxDegree=2) where {T,N,D}

    function flagPow(F, i)
        if i == 0
            return 1 // 1 * BinaryTreeFlag()
        elseif i == 1
            return 1 // 1 * F
        else
            return (1 // 1 * F) * flagPow(F, i - 1)
        end
    end

    fVars = Any[flagPow(G1, i) for i = 0:maxDegree]


    quadModuleModel = RazborovModel{T,N,D}(m)
    computeRazborovBasis!(quadModuleModel, razLevel - size(G1))

    function approxStep(minEdge, maxEdge, lowerBound=true)
        @show minEdge, maxEdge
        m.objective = nothing


        lm = addInequality!(m, G1 - minEdge * one(T), deepcopy(quadModuleModel))
        rm = addInequality!(m, maxEdge * one(T) - G1, deepcopy(quadModuleModel))

        computeSDP!.(m.subModels[end-1:end])


        jM = buildJuMPModel(m)
        set_optimizer(jM.model, solver)

        fact = lowerBound ? -1 : 1
        goalCoeff = Dict()
        for (A, d) in G2.coeff
            goalCoeff[A.tree] = get(goalCoeff, A.tree, 0) - fact * d
        end

        degVars = [1 * @variable(jM.model, base_name = "deg$(i-1)") for i in eachindex(fVars)]

        intF = AffExpr()
        for i in eachindex(fVars)
            degI = degVars[i]
            push!(degVars, degI)
            f = degI * fact
            for (A, d) in fVars[i].coeff
                goalCoeff[A.tree] = get(goalCoeff, A.tree, 0) + d * f
            end
            intF += fact * degI * (Float64(maxEdge)^i - Float64(minEdge)^i) / i
        end

        for (G, c) in jM.variables
            @constraint(jM.model, c == get(goalCoeff, G.tree, 0))
        end

        objective = fact * intF
        if lowerBound
            @objective(jM.model, Max, objective)
        else
            @objective(jM.model, Min, objective)
        end

        @time optimize!(jM.model)

        deleteat!(m.subModels, length(m.subModels))
        deleteat!(m.subModels, length(m.subModels))


        status = termination_status(jM.model)

        if status in [MOI.OPTIMAL, MOI.SLOW_PROGRESS, MOI.ALMOST_OPTIMAL]
            status == MOI.SLOW_PROGRESS && @warn "Slow progress!"
            coeff = deepcopy([value(degVars[i+1]) for i = 0:length(fVars)-1])

            lastPos = (maxEdge, sum(coeff[i+1] * maxEdge^i for i = 0:length(fVars)-1))

            firstPos = (minEdge, sum(coeff[i+1] * minEdge^i for i = 0:length(fVars)-1))

            return [firstPos, lastPos]
        else
            lastPos = nothing
            @error "Could not solve, $status"
            return [(minEdge, fact), (maxEdge, fact)]
        end
    end

    if lb
        lastPos = nothing
        lowerBounds = [approxStep(steps[i], steps[i+1], true) for i in 1:length(steps)-1]
    else
        lowerBound = 0
    end
    if ub
        lastPos = nothing
        upperBounds = [approxStep(steps[i], steps[i+1], false) for i in 1:length(steps)-1]
    else
        upperBound = 1
    end
    return lowerBounds, upperBounds
end

## Multiple levels of profile
profiles = Dict()

if isfile("Certificates/profiles.jld2")
    @load "Certificates/profiles.jld2" profiles
end
##
CairoMakie.activate!(type="svg")

trees = generateAll(BinaryTreeFlag, 7)
trees = labelCanonically.(trees)

interestingProfiles = [
    (trees[5], trees[14]),
    (trees[7], trees[14]),
    (trees[10], trees[14]),
    (trees[16], trees[14]),
    (trees[6], trees[9]),
    (trees[9], trees[13]),
    (trees[9], trees[12]),
    (trees[12], trees[13]),
    (trees[5], trees[15]),
]

for (t1, t2) in interestingProfiles

    # fig = Figure()
    # ax = Axis(fig[1, 1], xlabel="$t1", ylabel="$t2")#,aspect = DataAspect(),  limits=(0,1,0,1))


    if !haskey(profiles, (t1, t2))
        profiles[(t1, t2)] = Dict()
    end

    for lvl in max(size(t1), size(t2)):maxLvl
        if haskey(profiles[(t1, t2)], lvl)
            continue
        end


        if isfile("Models/modelLvl$lvl.jld2")
            @load "Models/modelLvl$lvl.jld2" m
        else
            rM = addRazborovBlock!(m, lvl)
            computeSDP!(m)

            !isdir("Models") && mkdir("Models")
            @save "Models/modelLvl$lvl.jld2" {compress = true} m
        end

        @show modelSize(m)
        T1 = 1 // 1 * t1
        T2 = 1 // 1 * t2

        m.objective = 1 * t1
        jm = buildJuMPModel(m)
        set_optimizer(jm.model, solver)
        optimize!(jm.model)
        @show termination_status(jm.model)
        t1Min = max(0 // 1, -rationalize(objective_value(jm.model) + 2e-5; tol=1e-5))

        m.objective = -1 * t1
        jm = buildJuMPModel(m)
        set_optimizer(jm.model, solver)
        optimize!(jm.model)
        @show termination_status(jm.model)
        t1Max = min(1 // 1, rationalize(objective_value(jm.model) + 2e-5; tol=1e-5))

        @time lower, upper = boundProfileQuantumFlag(m, T1, T2, lvl; maxDegree=1, steps=range(t1Min, t1Max, 100), solver=solver)

        profiles[(t1, t2)][lvl] = (lower, upper)

        # xPos = vcat([[l[1][1], l[2][1]] for l in lower]...)
        # yLower = vcat([[l[1][2], l[2][2]] for l in lower]...)
        # yUpper = vcat([[u[1][2], u[2][2]] for u in upper]...)
        # lines!(xPos, yUpper; color=:black)
        # band!(xPos, yLower, yUpper; color=(:black, 0.2))
        # lines!(xPos, yLower; color=:black)
        # lines!([xPos[1], xPos[1]], [yUpper[1], yLower[1]], color=:black)
        # lines!([xPos[end], xPos[end]], [yUpper[end], yLower[end]], color=:black)
        # display(fig)
        @save "Certificates/profiles.jld2" profiles
    end


end