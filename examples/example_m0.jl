using WignerFamilies
using PyPlot 
using OffsetArrays


j₂ = big(1000)
j₃ = big(1000)
m₁ = big(0)
m₂ = big(0)
m₃ = big(0)
js, w3j = WignerFamilies.classical_wigner3j_m0(BigFloat, j₂, j₃, m₂, m₃)


##
plt.clf()
plot(js, Float64.(w3j))
using WignerSymbols
plt.plot(js, 
    [wigner3j(Float64, n, j₂, j₃, 0, 0) for n in js], alpha=0.5)
plt.gcf()

##

plt.clf()
js = collect(w.nₘᵢₙ:w.nₘₐₓ)
using WignerSymbols
exact = OffsetArray([wigner3j(Float64, n, w.j₂, w.j₃, 0, 0) for n in js], w.nₘᵢₙ:w.nₘₐₓ)

plt.plot(Float64.(exact .- w3j))
# plt.plot(js, 
#     , alpha=0.5)
plt.gcf()

using Test
@test isapprox(exact, w3j)

