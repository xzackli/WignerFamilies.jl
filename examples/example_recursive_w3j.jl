using PyPlot 
using WignerFamilies

# compute all j₂=100, j₃=60, m₂=70, m₃=-55, m₁=-m₂-m₃
j₂, j₃, m₂, m₃ = 100, 60, 70, -55
m₁ = -m₂ - m₃
j_array, w3j = nonclassical_wigner3j(j₂, j₃, m₂, m₃)

plt.clf()
plt.plot(j_array, w3j .* 1000)
plt.xlabel(raw"$j$")
plt.ylabel(raw"$f(j) \, \times \, 10^3$")
plt.ylim(-20,20)
plt.axvline(49, ls="--", lw=1)
plt.axvline(98, ls="--", lw=1)
plt.tight_layout()
# plt.savefig("examples/luscombe_and_luban_1998.svg")
gcf()
