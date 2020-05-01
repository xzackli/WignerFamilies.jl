using WignerFamilies
using PyPlot 
using OffsetArrays

w = WignerF(300, 300, 0, 0)
start = Int( (w.nₘᵢₙ + w.nₘₐₓ) / 2 )
start = iseven(start) ? start : start + 1
println(w.nₘᵢₙ, " ", start, " ", w.nₘₐₓ)
w3j = OffsetArray(zeros(Float64, length(w.nₘᵢₙ:w.nₘₐₓ)), Int(w.nₘᵢₙ):Int(w.nₘₐₓ))
WignerFamilies.f_to_min_m0!(w, start, w3j)
WignerFamilies.f_to_max_m0!(w, start, w3j)

norm = WignerFamilies.normalization(w, w3j)
w3j ./= norm
# fjmax_sgn = iseven(w.j₂ + w.j₃ + w.m₂ + w.m₃) ? 1 : -1
if sign(w3j[w.nₘₐₓ]) != WignerFamilies.f_jmax_sgn(w)
    w3j .*= -1
end

##
using WignerSymbols
plt.clf()
js = collect(w.nₘᵢₙ:w.nₘₐₓ)
plot(js, Float64.(w3j))
exact = [wigner3j(BigFloat, n, w.j₂, w.j₃, 0, 0) for n in js]
plt.plot(js, Float64.(exact), alpha=0.5)
# plt.xlim(540, 560)
# plt.ylim(-0.002, 0.002)
plt.gcf()

##

plt.clf()
js = collect(w.nₘᵢₙ:w.nₘₐₓ)
using WignerSymbols
exact = OffsetArray([wigner3j(BigFloat, n, w.j₂, w.j₃, 0, 0) for n in js], w.nₘᵢₙ:w.nₘₐₓ)

plt.plot(Float64.(exact .- w3j))
# plt.plot(js, 
#     , alpha=0.5)
plt.gcf()