using WignerFamilies
using PyPlot 
using WignerSymbols

j₂, j₃, m₂, m₃ = 0, 10, 0, 0
w = WignerF(Float64, j₂, j₃, m₂, m₃)
println((w.nₘᵢₙ, w.nₘₐₓ))
w3j = get_wigner_array(w)
WignerFamilies.wigner3j_f!(w, w3j)
j_array = collect(eachindex(w3j))
w3j

##

plt.clf()
plot(j_array, Float64.(w3j.symbols))
plt.plot(j_array, 
    [Float64(wigner3j(BigFloat, n, j₂, j₃, 0, 0)) for n in j_array], alpha=0.5)
plt.gcf()



##


for j₂ in 1:100
    j₃, m₂, m₃ = 10, 0, 0
    w3j = wigner3j_f(Float64, j₂, j₃, m₂, m₃)
    j_array = collect(eachindex(w3j))
    for j in j_array
        @assert w3j[j] ≈ Float64(WignerSymbols.wigner3j(BigFloat, j, j₂, j₃, -m₂-m₃, m₂))
    end
end

