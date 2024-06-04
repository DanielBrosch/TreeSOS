using FlagSOS, JLD2, JuMP, LinearAlgebra, ProgressMeter, CairoMakie

lvl = 8

function latexify(T)
    tmp = "$(FlagSOS.labelLex(T.tree))"
    tmp = replace(tmp, "•" => "\\bullet")
    return L"%$tmp"
end

trees = generateAll(BinaryTreeFlag, lvl)
trees = labelCanonically.(trees)
t1 = trees[5]
t2 = trees[14]
fig = Figure()
ax = Axis(fig[1, 1], xlabel=latexify(t1), ylabel=latexify(t2))

profiles = Dict()
if isfile("Certificates/profiles.jld2")
    @load "Certificates/profiles.jld2" profiles
end

profs = collect(profiles[(t1, t2)])
sort!(profs, by = x -> x[1])

for (i, (lower, upper)) in profs
    i > 8 && continue
    xPos = vcat([[l[1][1], l[2][1]] for l in lower]...)
    yLower = vcat([[l[1][2], l[2][2]] for l in lower]...)
    yUpper = vcat([[u[1][2], u[2][2]] for u in upper]...)
    lines!(xPos, yUpper; color=:black, linewidth=0.5)
    band!(xPos, yLower, yUpper; color=(:black, 0.2))
    lines!(xPos, yLower; color=:black, linewidth=0.5)
    lines!([xPos[1], xPos[1]], [yUpper[1], yLower[1]], color=:black, linewidth=0.5)
    lines!([xPos[end], xPos[end]], [yUpper[end], yLower[end]], color=:black, linewidth=0.5)
end


cD = t -> t^4 + (1 - t)^4 + 4 * (t * (1 - t)^3 + t^3 * (1 - t))
bD = t -> 20 * t^3 * (1 - t)^3
steps = range(0,1,500)
lines!(cD.(steps), bD.(steps))

t = 1 // 4

m = FlagModel{BinaryTreeFlag,:limit,Rational{Int}}()
rM = addRazborovBlock!(m, lvl)


qM = FlagSOS.EqualityModule{BinaryTreeFlag,BinaryTreeFlag,:limit,Rational{Int}}(1 // 1 * t1 - cD(t) * one(t1))

Fs = generateAll(BinaryTreeFlag, lvl - size(t1))[4:5]

for F in Fs
    push!(qM.basis, F)
end

push!(m.subModels, qM)

computeSDP!(m)

m.objective = -1 // 1 * t2
type = Float64

jm = buildJuMPModel(m, Dict(), GenericModel{type}(), false)
set_optimizer(jm.model, solver)

optimize!(jm.model)

## Rigorous solution

using LinearAlgebra

ex = FlagSOS.roundResults(m, jm...; prec=1 // 10^8)

!isdir("Certificates") && mkdir("Certificates")
open("Certificates/NonconvexityCat4E6.txt", "w") do io
    res = FlagSOS.verifySOS(m, ex; io=io)
    @show roundedUpperBound = res.coeff[one(BinaryTreeFlag)]
    @show floatBound = Float64(roundedUpperBound, RoundUp)
end


rig = scatter!([cD(t)], [objective_value(jm.model)], label = "rigorous bound", color=:red)
Legend(fig[1,2],
    [rig],
    ["rigorous upper bound"])

!isdir("Pictures") && mkdir("Pictures")
save("Pictures/nonconvexity.svg", fig)
save("Pictures/nonconvexity.pdf", fig)