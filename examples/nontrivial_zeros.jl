using PyPlot 
using WignerFamilies
using HalfIntegers
# compute all j₂=100, j₃=60, m₂=70, m₃=-55, m₁=-m₂-m₃
# bottom of page 14 of Raynal et al.
# On the zeros of 3j coefficients: polynomial degree versus recurrence order

X₀, Y₀ = 43, 25
# X₀, Y₀ = 109, 63
X(n) = iszero(n) ? X₀ : 7X(n-1) + 12Y(n-1)
Y(n) = iszero(n) ? Y₀ : 4(n-1) + 7Y(n-1)

n_seq = 1
# j₂, j₃, m₂, m₃ = 100, 60, 70, -55
# j₂, j₃, m₂, m₃ = (2X(n_seq)+1)/6, (X(n_seq)+2)/6, HalfInt(3/2), HalfInt(3/2)
# j₂, j₃, m₂, m₃ = 4, 4, 2, -2

j₂, j₃, m₂, m₃ = HalfInt(23/2), HalfInt(21/2), HalfInt(3/2), HalfInt(3/2)

m₁ = -m₂ - m₃

w = WignerF(BigFloat, j₂, j₃, m₂, m₃)


nmid = Int(ceil((w.nₘᵢₙ + w.nₘₐₓ)/2))
nmid = iseven(nmid) ? nmid : nmid + 1
println(w.nₘᵢₙ, " ", nmid, " ", w.nₘₐₓ)
w3j = get_wigner_array(w)

n₊ = WignerFamilies.rψ!(w, nmid, w3j)
n₋ = WignerFamilies.sψ!(w, nmid, w3j)

for k in (n₊+1):w.nₘₐₓ  # product the ratios upwards
    w3j[k] *= w3j[k-1]
end
for k in (n₋-1):-1:w.nₘᵢₙ  # product the ratios downwards
    w3j[k] *= w3j[k+1]
end

println(n₋, ", ", n₊)

ψn₋ = w3j[n₋]
ψn₊ = w3j[n₊]

WignerFamilies.ψauxplus!(w, n₋, n₊, w3j)

if n₋ > w.nₘᵢₙ
    scale_middle = ψn₋ / w3j[n₋]
    w3j[w.nₘᵢₙ:(n₋+1)] .*= scale_middle
end
if n₊ < w.nₘₐₓ
    scale_end = w3j[n₊] / ψn₊
    w3j[(n₊+1):w.nₘₐₓ] .*= scale_end
end

norm = WignerFamilies.normalization(w, w3j)
w3j ./= norm
if sign(w3j[w.nₘₐₓ]) != WignerFamilies.f_jmax_sgn(w)
    w3j .*= -1
end

##

plt.clf()
using WignerSymbols
exact = [wigner3j(BigFloat, n, w.j₂, w.j₃, m₁, m₂) for n in w.nₘᵢₙ:w.nₘₐₓ]
plt.plot(w.nₘᵢₙ:w.nₘₐₓ, Float64.(w3j))

plt.plot(w.nₘᵢₙ:w.nₘₐₓ, Float64.(exact), alpha=0.5)
plt.xlabel(raw"$j$")
plt.ylabel(raw"$f(j) \, \times \, 10^3$")
plt.tight_layout()
# plt.xlim(100,120)
gcf()


##
plt.clf()
plt.plot( Float64.(exact ./ w3j.parent))
plt.gcf()

##