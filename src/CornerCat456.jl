using FlagSOS, JLD2, JuMP, LinearAlgebra, ProgressMeter, CairoMakie

lvl = 6

trees = generateAll(BinaryTreeFlag, lvl)
trees = labelCanonically.(trees)
E6 = trees[14]

for Cat in trees[[5, 7, 10]]

    # Corner: (size(Cat)+1)//2^(size(Cat)-1), 5//16

    m = FlagModel{BinaryTreeFlag,:limit,Rational{Int}}()
    rM = addRazborovBlock!(m, lvl)

    blocks = sort!(collect(keys(rM.basis)), by=x -> string(x), rev=true)
    deleteBlocks = [1, 2, 3, 4, 7, 8, 9]

    for b in deleteBlocks
        delete!(rM.basis, blocks[b])
    end

    qM = FlagSOS.EqualityModule{BinaryTreeFlag,BinaryTreeFlag,:limit,Rational{Int}}(1 // 1 * Cat - (size(Cat) + 1) // 2^(size(Cat) - 1) * one(Cat))

    Fs = generateAll(BinaryTreeFlag, lvl - size(Cat))
    push!(qM.basis, one(BinaryTreeFlag))
    push!(m.subModels, qM)

    computeSDP!(m)


    m.objective = -30 * E6 # scaled such that all coefficients are integers
    type = Float64

    jm = buildJuMPModel(m, Dict(), GenericModel{type}(), false)

    set_optimizer(jm.model, solver)

    optimize!(jm.model)
    @show termination_status(jm.model)
    @show objective_value(jm.model)

    # Exact solution

    ex = [Dict(), Dict()]

    for (i, B) in enumerate(jm.blocks)
        for (mu, b) in B
            ex[i][mu] = round.(Int, value(b))
        end
    end


    !isdir("Certificates") && mkdir("Certificates")
    open("Certificates/MiddleCornerCat$(size(Cat)).txt", "w") do io
        res = FlagSOS.verifySOS(m, ex; io=io)
        res = 1 // res.coeff[E6] * res
        @show res
    end

end