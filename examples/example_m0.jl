using WignerFamilies
using PyPlot 
using OffsetArrays

j₂, j₃, m₂, m₃ = 100, 60, 0, 0
w = WignerF(Float64, j₂, j₃, m₂, m₃)
w3j = get_wigner_array(w)
@allocated WignerFamilies.classical_wigner3j_m0!(w, w3j)

##
plt.clf()
plot(js, Float64.(w3j.symbols))
using WignerSymbols
plt.plot(js, 
    [Float64(wigner3j(BigFloat, n, j₂, j₃, 0, 0)) for n in js], alpha=0.5)
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

