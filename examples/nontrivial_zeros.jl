using PyPlot 
using WignerFamilies

# compute all j₂=100, j₃=60, m₂=70, m₃=-55, m₁=-m₂-m₃
# bottom of page 14 of Raynal et al.
# On the zeros of 3j coefficients: polynomial degree versus recurrence order
X = 43
j₂, j₃, m₂, m₃ = (2X+1)/6, (X+2)/6, HalfInt(3/2), HalfInt(3/2)
m₁ = -m₂ - m₃
w = WignerF(Float64, j₂, j₃, m₂, m₃)

start = Int( ceil((w.nₘᵢₙ + w.nₘₐₓ) / 2) )
start = iseven(start) ? start : start + 1
println(w.nₘᵢₙ, " ", start, " ", w.nₘₐₓ)
w3j = get_wigner_array(w)
w3j .= NaN

println(WignerFamilies.rψ!(w, 20, w3j))
# WignerFamilies.f_to_min_m0!(w, start, w3j)
# WignerFamilies.f_to_max_m0!(w, start, w3j)

print(w3j)


##


plt.clf()
using WignerSymbols
exact = [wigner3j(BigFloat, n, w.j₂, w.j₃, m₁, m₃) for n in w.nₘᵢₙ:w.nₘₐₓ]
println( count(iszero, exact), " ", w.nₘᵢₙ:w.nₘₐₓ)
plt.plot(Float64.(exact))
plt.gcf()


##
plt.clf()
plt.plot(j_array, w3j .* 1000)
plt.xlabel(raw"$j$")
plt.ylabel(raw"$f(j) \, \times \, 10^3$")
# plt.ylim(-20,20)
# plt.axvline(49, ls="--", lw=1)
# plt.axvline(98, ls="--", lw=1)
plt.tight_layout()
# plt.savefig("examples/luscombe_and_luban_1998.svg")
gcf()
