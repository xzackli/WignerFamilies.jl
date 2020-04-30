using PyPlot 
using WignerFamilies

js, w3j = nonclassical_wigner3j(100, 60, 70, -55)

plt.clf()
plt.plot(js, w3j .* 1000)
plt.xlabel(raw"$j$")
plt.ylabel(raw"$f(j) \, \times \, 10^3$")
# plt.yticks([-20, -15, -10, -5, 0, 5, 10, 15, 20])
plt.ylim(-20,20)
plt.axvline(49, ls="--", lw=1)
plt.axvline(98, ls="--", lw=1)
plt.tight_layout()
plt.savefig("examples/luscombe_and_luban_1998.png")
gcf()