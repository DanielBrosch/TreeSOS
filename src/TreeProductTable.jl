using FlagSOS, JLD2, JuMP, LinearAlgebra, ProgressMeter, CairoMakie

lvl = 6

trees = generateAll(BinaryTreeFlag, lvl)
trees = labelCanonically.(trees)

m = FlagModel{BinaryTreeFlag,:limit,Rational{Int}}()
rM = addRazborovBlock!(m, lvl)
rM2 = addRazborovBlock!(m, lvl - 1)

modelSize(rM)
modelSize(rM2)

modelSize(m)
computeSDP!(m)

##
allProducts = []


blocks = vcat(collect(rM.basis), collect(rM2.basis))
sort!(blocks, by=x -> size(x[1]))

!isdir("GeneratedLatex") && mkdir("GeneratedLatex")
open("GeneratedLatex/treeProducts.tex", "w") do io
    println(io, "{\\allowdisplaybreaks\\tiny")
    for (T, b) in blocks
        bUnsorted = deepcopy(b)
        sort!(b, by = x->(size(x), string(x)))
        @show b
        razM = haskey(rM.basis, T) ? rM : rM2
        println(io, "Type \\resizebox{.06\\textwidth}{!}{" * FlagSOS.printTreeLatex(T) * "}")
        println(io, "\\begin{align*}")
        coveredInds = Int[]
        for (it, b1) in enumerate(b)
            for (jt, b2) in enumerate(b)
                # i > j && continue
                i = findfirst(x->x==b1, bUnsorted)
                j = findfirst(x->x==b2, bUnsorted)
                # if CartesianIndex(j, i) != findfirst(x -> x == razM.blockSymmetry[T].pattern[i, j], razM.blockSymmetry[T].pattern)
                #     continue
                # end
                p = razM.blockSymmetry[T].pattern[i, j]
                if in(p, coveredInds)
                    continue 
                end
                push!(coveredInds, p)
                # print(io, "\$\$" * FlagSOS.printTreeLatex(b1) * " \\cdot " * FlagSOS.printTreeLatex(b2) * " = ")
                print(io, "\\llbracket $b1" * " \\cdot " * "$b2\\rrbracket" * " &= ")
                first = true
                count = 0
                tmp = prod(razM.sdpData) do (G, B)
                    if haskey(B, T)
                        c = B[T][i, j]
                        iszero(c) && return ""
                        fraction = ""
                        if c.den == 1
                            if c.num != 1
                                fraction = "$(c.num)"
                            end
                        else
                            fraction = "\\textstyle\\tfrac{$(c.num)}{$(c.den)}"
                        end
                        if !first
                            fraction = "+ " * fraction
                        end
                        if count == 2
                            fraction = "\\\\&\\qquad " *fraction
                            count -= 2
                        end
                        first = false
                        count += 1
                        return fraction * "$G"
                        # return fraction * FlagSOS.printTreeLatex(G)
                    else
                        return ""
                    end
                end
                println(io, tmp * "\\\\")
            end
        end

        println(io, "\\end{align*}")
    end

    println(io, "Quotient:")
    println(io, "\\begin{align*}")
    for F in rM.quotient
        @show F
        neg = findfirst(x -> x == -1, F.coeff)
        FPos = F + neg
        @show FPos
        first = true
        tmp = prod(FPos.coeff) do (G, c)
            iszero(c) && return ""
            fraction = ""
            if c.den == 1
                if c.num != 1
                    fraction = "$(c.num)"
                end
            else
                fraction = "\\tfrac{$(c.num)}{$(c.den)}"
            end
            if !first
                fraction = "+ " * fraction
            end
            first = false
            return fraction * "$G"
        end
        println(io, "$neg &= \\bullet \\cdot $neg \\\\&= " * tmp * "\\\\")
    end
    println(io, "\\end{align*}")

    println(io, "}")
end