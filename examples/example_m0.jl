using WignerFamilies
using PyPlot 

j₂, j₃, m₂, m₃ = 1, 0, 0, 0
w = WignerF(Float64, j₂, j₃, m₂, m₃)
println((w.nₘᵢₙ, w.nₘₐₓ))
w3j = get_wigner_array(w)
WignerFamilies.wigner3j_f!(w, w3j)
# j_array = collect(eachindex(w3j))
w3j


##

wigner3j_f(j₂, j₃, 0, 0)

##
plt.clf()
plot(j_array, Float64.(w3j.symbols))
using WignerSymbols
plt.plot(j_array, 
    [Float64(wigner3j(BigFloat, n, j₂, j₃, 0, 0)) for n in j_array], alpha=0.5)
plt.gcf()

##

plt.clf()
js = collect(w.nₘᵢₙ:w.nₘₐₓ)
using WignerSymbols
exact = [wigner3j(Float64, n, w.j₂, w.j₃, 0, 0) for n in js]
plt.plot(Float64.(exact .- w3j.symbols))
plt.gcf()

using Test
@test isapprox(exact, w3j.symbols)

