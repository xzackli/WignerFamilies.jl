# WignerFamilies

[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://xzackli.github.io/WignerFamilies.jl/dev)
[![Build Status](https://github.com/xzackli/WignerFamilies.jl/workflows/CI/badge.svg)](https://github.com/xzackli/WignerFamilies.jl/actions)
[![Coverage](https://codecov.io/gh/xzackli/WignerFamilies.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/xzackli/WignerFamilies.jl)
<!-- ![Lifecycle](https://img.shields.io/badge/lifecycle-experimental-orange.svg) -->
<!-- [![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://xzackli.github.io/WignerFamilies.jl/stable) -->

This package implements methods described in [Luscombe and Luban 1998](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.57.7274), based on the work of [Schulten and Gordon 1961](https://aip.scitation.org/doi/10.1063/1.522426), for generating families of Wigner 3j and 6j symbols by recurrence relation. These exact methods are orders of magnitude more efficient than than strategies like [prime factorization](https://github.com/Jutho/WignerSymbols.jl) for problems which require every non-trivial symbol in a family, and really shine for large quantum numbers. WignerFamilies.jl is thread-safe and **very fast**, beating the standard Fortran routine DRC3JJ from SLATEC by a factor of 2-4 (see [notebook](https://nbviewer.jupyter.org/github/xzackli/WignerFamilies.jl/blob/master/test/benchmarks.ipynb)).

## Installation

```julia
using Pkg
Pkg.add("WignerFamilies")
```

## Usage
WignerFamilies.jl currently computes the nontrivial 3j symbols over `j` with the other 
quantum numbers fixed, in the family of symbols,

<p align="center">
<img width=40% src="docs/src/assets/fsymbol.png">
</p>

It exposes `wigner3j_f(j₂, j₃, m₂, m₃)` which returns a simple wrapper around a vector of 
the type`WignerSymbolVector`. This vector contains the computed symbols, indexed by the 
quantum number `j`. The type supports 
[half-integer](https://github.com/sostock/HalfIntegers.jl) quantum numbers as indices.

```julia
using WignerFamilies

# wigner3j for all j fixing j₂=100, j₃=60, m₂=70, m₃=-55, m₁=-m₂-m₃
w3j = wigner3j_f(100, 60, 70, -55) 
js = collect(eachindex(w3j))  # indices are the quantum numbers
plot(js, w3j.symbols)   # you can get the underlying array with w3j.symbols
```

This generates the symbols in Figure 1 of Luscombe and Luban 1998.
<p align="center">
<img width=90% src="docs/src/assets/luscombe_and_luban_1998.svg">
</p>

## Advanced Use

One can compute symbols in a fully non-allocating way by using the mutating `wigner3j_f!`. This requires 
one to initialize a `WignerF` struct describing the problem, put a wrapper around the piece of memory
you want to deposit the symbols using WignerSymbolVector, and then calling the mutating method.

```julia
j₂, j₃, m₂, m₃ = 100, 100, -2, 2
w = WignerF(Float64, j₂, j₃, m₂, m₃)
buffer = zeros(Float64, length(w.nₘᵢₙ:w.nₘₐₓ))
w3j = WignerSymbolVector(buffer, w.nₘᵢₙ:w.nₘₐₓ)
WignerFamilies.wigner3j_f!(w, w3j)
```

This is the best way to use to get symbols if you're using them in a tight loop, since allocations really hurt in those cases.

```julia
using BenchmarkTools
@btime WignerFamilies.wigner3j_f!(w, w3j)
```
```
1.660 μs (0 allocations: 0 bytes)
```
