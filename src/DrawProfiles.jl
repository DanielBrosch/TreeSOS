using FlagSOS, JLD2, JuMP, LinearAlgebra, ProgressMeter, CairoMakie

function latexify(T)
    tmp = "$(FlagSOS.labelLex(T.tree))"
    tmp = replace(tmp, "•" => "\\bullet")
    return L"%$tmp"
end


@load "Certificates/profiles.jld2" profiles

trees = generateAll(BinaryTreeFlag, 7)
trees = labelCanonically.(trees)

profileFigs = Dict()

for (t1, t2) in keys(profiles)

    fig = Figure()
    # ax = Axis(fig[1, 1], xlabel="$(FlagSOS.labelLex(t1.tree))", ylabel="$(FlagSOS.labelLex(t2.tree))")#,aspect = DataAspect(),  limits=(0,1,0,1))
    ax = Axis(fig[1, 1], xlabel=latexify(t1), ylabel=latexify(t2))#,aspect = DataAspect(),  limits=(0,1,0,1))

    profs = collect(profiles[(t1, t2)])
    sort!(profs, by=x -> x[1])

    for (i, (lower, upper)) in profs
        xPos = vcat([[l[1][1], l[2][1]] for l in lower]...)
        yLower = vcat([[l[1][2], l[2][2]] for l in lower]...)
        yUpper = vcat([[u[1][2], u[2][2]] for u in upper]...)
        lines!(xPos, yUpper; color=:black, linewidth=0.5)
        band!(xPos, yLower, yUpper; color=(:black, 0.2))
        lines!(xPos, yLower; color=:black, linewidth=0.5)
        lines!([xPos[1], xPos[1]], [yUpper[1], yLower[1]], color=:black, linewidth=0.5)
        lines!([xPos[end], xPos[end]], [yUpper[end], yLower[end]], color=:black, linewidth=0.5)
        # display(fig)
    end
    t1n = findfirst(x -> x == t1, trees)
    t2n = findfirst(x -> x == t2, trees)
    !isdir("Pictures") && mkdir("Pictures")
    save("Pictures/profile$(t1n)_$t2n.svg", fig)
    save("Pictures/profile$(t1n)_$t2n.pdf", fig)
    profileFigs[(t1, t2)] = fig
end

## Profiles with lower bounds 

fig = Figure()
trees = generateAll(BinaryTreeFlag, 7)
trees = labelCanonically.(trees)

lowerFigs = Dict()

t2 = trees[14]
t1s = trees[[5, 7, 10]]

for (k, t1) in enumerate(t1s)

    if k == 1
        ax = Axis(fig[1, k], ylabel=latexify(t2))
    else
        ax = Axis(fig[1, k])
    end

    profs = collect(profiles[(t1, t2)])
    sort!(profs, by=x -> x[1])

    for (i, (lower, upper)) in profs
        xPos = vcat([[l[1][1], l[2][1]] for l in lower]...)
        yLower = vcat([[l[1][2], l[2][2]] for l in lower]...)
        yUpper = vcat([[u[1][2], u[2][2]] for u in upper]...)
        lines!(ax, xPos, yUpper; color=:black, linewidth=0.5)
        band!(ax, xPos, yLower, yUpper; color=(:black, 0.2))
        lines!(ax, xPos, yLower; color=:black, linewidth=0.5)
        lines!(ax, [xPos[1], xPos[1]], [yUpper[1], yLower[1]], color=:black, linewidth=0.5)
        lines!(ax, [xPos[end], xPos[end]], [yUpper[end], yLower[end]], color=:black, linewidth=0.5)
        # display(fig)
        t1n = findfirst(x -> x == t1, trees)
        t2n = findfirst(x -> x == t2, trees)
    end


    ax = Axis(fig[2, k], xlabel=latexify(t1))
    sharpInd = [7, 8, 8][k]
    @show (t1, t2)
    @show sharpInd 
    @show keys(profiles[(t1,t2)])
    firstLower = profiles[(t1, t2)][min(sharpInd, maximum(keys(profiles[(t1,t2)])))][1]
    firstYLower = vcat([[l[1][2], l[2][2]] for l in firstLower]...)

    profs = collect(profiles[(t1, t2)])
    sort!(profs, by=x -> x[1])

    for (i, (lower, upper)) in profs
        xPos = vcat([[l[1][1], l[2][1]] for l in lower]...)
        yLower = vcat([[l[1][2], l[2][2]] for l in lower]...) - firstYLower
        yUpper = vcat([[u[1][2], u[2][2]] for u in upper]...)

        lines!(ax, xPos, yLower; label="$i", color=:black, linewidth=1)
    end
    # display(fig)
    lowerFigs[(t1, t2)] = fig
end


colsize!(fig.layout, 1, Aspect(1, 1.0))
colsize!(fig.layout, 2, Aspect(1, 1.0))
colsize!(fig.layout, 3, Aspect(1, 1.0))
Label(fig[1, 3, Right()], "Profiles", rotation=-pi / 2, fontsize=15, padding=(15.0f0, 0.0f0, 0.0f0, 0.0f0))
Label(fig[2, 3, Right()], "Lower bound differences", rotation=-pi / 2, fontsize=15, padding=(15.0f0, 0.0f0, 0.0f0, 0.0f0))
resize_to_layout!(fig)
# display(fig)


!isdir("Pictures") && mkdir("Pictures")
save("Pictures/caterpillars.svg", fig)
save("Pictures/caterpillars.pdf", fig)

## Caterpillar 4 vs balanced

trees = generateAll(BinaryTreeFlag, 7)
trees = labelCanonically.(trees)

t1 = trees[5]
t2 = trees[14]

fig = Figure()


xL = L"C_4 = (\bullet(\bullet(\bullet\bullet)))"
yL = L"E_6 = ((\bullet(\bullet\bullet))(\bullet(\bullet\bullet)))"


ax = Axis(fig[1, 1], xlabel=xL, ylabel=yL, xautolimitmargin=(0.05f0, 0.09f0), yautolimitmargin=(0.05f0, 0.12f0))

firstUpper = profiles[(t1, t2)][min(8, maximum(keys(profiles[(t1,t2)])))][2]
firstYUpper = vcat([[l[1][2], l[2][2]] for l in firstUpper]...)

profs = collect(profiles[(t1, t2)])
sort!(profs, by=x -> x[1])

profs = profs[end:end]

for (i, (lower, upper)) in profs
    xPos = vcat([[l[1][1], l[2][1]] for l in lower]...)
    yLower = vcat([[l[1][2], l[2][2]] for l in lower]...)
    yUpper = vcat([[u[1][2], u[2][2]] for u in upper]...)
    lines!(xPos, yUpper; color=:black, linewidth=0.5)
    band!(xPos, yLower, yUpper; color=(:black, 0.2))
    lines!(xPos, yLower; color=:black, linewidth=0.5)
    lines!([xPos[1], xPos[1]], [yUpper[1], yLower[1]], color=:black, linewidth=0.5)
    lines!([xPos[end], xPos[end]], [yUpper[end], yLower[end]], color=:black, linewidth=0.5)
end

# display(fig)
t1n = findfirst(x -> x == t1, trees)
t2n = findfirst(x -> x == t2, trees)
fy = t -> 20 * t^3 * (1 - t)^3
fx = t -> t^4 + (1 - t)^4 + 4 * (t * (1 - t)^3 + t^3 * (1 - t))

gy = t -> 10 / 31 * (t^6 + (1 - t)^6) + 20 * t^3 * (1 - t)^3
gx = t -> begin
    factorial(4) / 2 * prod((2^j - 1)^(-1) for j = 1:2) * *((t^4 + (1 - t)^4) / (2^3 - 1) + t^3 * (1 - t) + (1 - t)^3 * t)
end

tRange = range(1 // 2, 1 // 1, 500)
lines!(fx.(tRange), fy.(tRange), color=:red)


text!(5 // 8 + 0.015, 5 // 16 - 0.01; text=L"\mathrm{DCat}^{1/2,1/2}_\infty:\left(\frac{5}{8},\enspace \frac{5}{16}\right)", color=:blue)
scatter!([5 // 8], [5 // 16], color=:blue, markersize=10)
text!(1 // 1 - 0.04, 0 // 1 + 0.015; text=L"C_{\infty}:\left(\frac{1}{1},\enspace \frac{0}{1}\right)", color=:blue)
scatter!([1 // 1], [0 // 1], color=:blue, markersize=10)
text!(4 // 7 - 0.01, 10 // 31 + 0.01; text=L"E_{\infty}:\left(\frac{4}{7},\enspace \frac{10}{31}\right)", color=:blue)
scatter!([4 // 7], [10 // 31], color=:blue, markersize=10)


text!(0.66, 0.27; text=L"\mathrm{DCat}^{t,1-t}:\left(20t^3(1-t)^3, t^4+(1-t^4)+4(t^3(t-1)+t(t-1)^3)\right)", color=:red, rotation=1.78 * pi)


!isdir("Pictures") && mkdir("Pictures")
save("Pictures/profile_marked$(t1n)_$t2n.svg", fig)
save("Pictures/profile_marked$(t1n)_$t2n.pdf", fig)
# fig
