using FlagSOS, JLD2, JuMP, LinearAlgebra, ProgressMeter, CairoMakie

!isdir("GeneratedLatex") && mkdir("GeneratedLatex")
open("GeneratedLatex/blockTable.tex", "w") do io
    println(io, "\\begin{tabular}{@{}ccc@{}}")
	println(io, "\\toprule")
    println(io, "Level & Block sizes & Sum \\\\\\midrule")
    for lvl in 4:maxLvl
        !isfile("Models/modelLvl$lvl.jld2") && continue
        print(io, "\$$lvl\$ & \$")

        @load "Models/modelLvl$lvl.jld2" m
        bs = modelSize(m).part
        for b in sort(unique(bs); rev=true)
            bc = count(x->x==b, bs)
            print(io, "$(b)_{$bc}")
        end
        println(io, "\$ & \$$(sum(bs))\$\\\\")
        # @show modelSize(m)
        # @show sum(modelSize(m))
    end
    println(io, "\\bottomrule")
    println(io, "\\end{tabular}")


end
