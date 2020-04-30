


##
using PyPlot 

w = WignerF(80, 10, -2, 2)
start = Int( (w.nₘᵢₙ + w.nₘₐₓ) / 2 )
start = iseven(start) ? start : start + 1
println(w.nₘᵢₙ, " ", start, " ", w.nₘₐₓ)
w3j = OffsetArray(zeros(BigFloat, length(w.nₘᵢₙ:w.nₘₐₓ)), w.nₘᵢₙ:w.nₘₐₓ)
f_to_min_m0!(w, start, w3j)
f_to_max_m0!(w, start, w3j)

norm = normalization(w, w3j)
w3j ./= norm
fjmax_sgn = iseven(w.j₂ + w.j₃ + w.m₂ + w.m₃) ? 1 : -1
if sign(w3j[w.nₘₐₓ]) != fjmax_sgn
    w3j .*= -1
end


plt.clf()
js = collect(w.nₘᵢₙ:w.nₘₐₓ)
plot(js, Float64.(w3j))
using WignerSymbols
plt.plot(js, 
    [wigner3j(Float64, n, w.j₂, w.j₃, 0, 0) for n in js], alpha=0.5)
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

