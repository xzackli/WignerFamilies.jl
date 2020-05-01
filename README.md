# WignerFamilies

[![Build Status](https://github.com/xzackli/WignerFamilies.jl/workflows/CI/badge.svg)](https://github.com/xzackli/WignerFamilies.jl/actions)
[![Coverage](https://codecov.io/gh/xzackli/WignerFamilies.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/xzackli/WignerFamilies.jl)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://xzackli.github.io/WignerFamilies.jl/dev)
![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg)
<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://xzackli.github.io/WignerFamilies.jl/stable) -->

This implements the methods described in [Luscombe and Luban 1998](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.57.7274), based on the work of [Schulten and Gordon 1961](https://aip.scitation.org/doi/10.1063/1.522426), for generating families of Wigner 3j and 6j symbols by recurrence relation. It also contains code implementing the magic square methods for Wigner 3j symbols in [Rasch and Yu 2012](https://epubs.siam.org/doi/abs/10.1137/S1064827503422932).

```julia
using WignerFamilies

# wigner3j for all j₁ fixing j₂=100, j₃=60, m₂=70, m₃=-55, m₁=-m₂-m₃
w3j = wigner3j_f(100, 60, 70, -55)  # outputs an OffsetArray for j₁'s nontrivial interval
j_array = collect(eachindex(w3j))
```

This generates Figure 1 of Luscombe and Luban 1998.
<p align="center">
![example plot](examples/luscombe_and_luban_1998.svg)
</p>
