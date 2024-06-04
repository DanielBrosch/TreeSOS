using FlagSOS, JLD2, JuMP, LinearAlgebra, ProgressMeter, CairoMakie

## Computes bounds on the inducibility of trees with up to maxLvl leaves, rounds them, and generates a table

for lvl in 4:maxLvl

    trees = generateAll(BinaryTreeFlag, lvl)[2:end]
    trees = labelCanonically.(trees)
    m = FlagModel{BinaryTreeFlag,:limit,Rational{Int}}()

    if isfile("Models/modelLvl$lvl.jld2")
        @load "Models/modelLvl$lvl.jld2" m
    else
        rM = addRazborovBlock!(m, lvl)
        computeSDP!(m)

        !isdir("Models") && mkdir("Models")
        @save "Models/modelLvl$lvl.jld2" {compress = true} m
    end

    modelSize(m) # lvl 9: 70₁54₃13₁₁9₁1₄₆

    #

    inducibilityBounds = Dict()
    inducibilityCertificates = Dict()

    irrationalTree = BinaryTree(BinaryTree(false), BinaryTree(BinaryTree(BinaryTree(false), BinaryTree(false)), BinaryTree(BinaryTree(false), BinaryTree(false))))

    @showprogress "Solving SDPs..." for T in trees
        m.objective = -1 * T
        jm = buildJuMPModel(m)

        set_optimizer(jm.model, solver)
        optimize!(jm.model)
        @show T
        # @show termination_status(jm.model)
        @show objective_value(jm.model)
        inducibilityBounds[T.tree] = objective_value(jm.model)
        inducibilityCertificates[T.tree] = jm
    end
    @show inducibilityBounds


    !isdir("Bounds") && mkdir("Bounds")
    @save "Bounds/InducibilityBounds$lvl.jld2" {compress = true} inducibilityBounds

    # Rounding results 

    roundedBounds = Dict()
    roundedCertificates = Dict()
    @showprogress "Rounding certificates..." for T in trees
        m.objective = -1 * T
        exT = FlagSOS.roundResults(m, inducibilityCertificates[T.tree]...; prec=1 // 10^8)
        roundedCertificates[T.tree] = exT
        bound = FlagSOS.verifySOS(m, exT; io=Base.DevNull())
        roundedBounds[T.tree] = bound
        @show inducibilityBounds[T.tree]
        display(convert(QuantumFlag{BinaryTreeFlag,Float64}, bound))

    end


    !isdir("Bounds") && mkdir("Bounds")
    !isdir("Certificates") && mkdir("Certificates")
    @save "Bounds/InducibilityBoundsRounded$lvl.jld2" {compress = true} roundedBounds
    @save "Certificates/InducibilityCertificatesRounded$lvl.jld2" {compress = true} roundedCertificates

end

## Inducibility table

@load "Bounds/InducibilityBoundsRounded$(maxLvl).jld2" roundedBounds
roundedBounds11 = roundedBounds
@load "Bounds/InducibilityBoundsRounded$(maxLvl-1).jld2" roundedBounds
roundedBounds10 = roundedBounds
# @load "Certificates/InducibilityBounds11.jld2" inducibilityBounds

e = one(BinaryTreeFlag)

# bounds = collect(inducibilityBounds)
# sort!(bounds, by=x -> (size(x[1]), string(FlagSOS.labelLex(x[1]))))


!isdir("GeneratedLatex") && mkdir("GeneratedLatex")
open("GeneratedLatex/InducibilityBounds.tex", "w") do io
    println(io, "{\\scriptsize\\allowdisplaybreaks\\begin{align*}")

    for (T, c) in roundedBounds11
        print(io, "$(FlagSOS.labelLex(T))&")
        tmp11 = 0.0
        # if haskey(roundedBounds11, T)
        c2 = Float64(roundedBounds11[T].coeff[e], RoundUp)
        tmp11 = round(c2, RoundUp; digits=7)
        print(io, "\\leq $(tmp11)")
        # else
        #     tmp11 = round(c, RoundUp; digits=7)
        #     print(io, "\\red{\\leq $(tmp11)}")
        # end
        if haskey(roundedBounds10, T)
            c2 = Float64(roundedBounds[T].coeff[e], RoundUp)
            tmp = round(c2, RoundUp; digits=7)
            if tmp11 <= tmp
                print(io, "\\leq ")
            else
                print(io, "\\red{\\leq} ")
            end
            print(io, "$(tmp)")
        end
        println(io, "\\\\")

    end
    println(io, "\\end{align*}}")
end