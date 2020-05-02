using WignerFamilies
using PyPlot 
using OffsetArrays
using BenchmarkTools

w = WignerF(Float64, 4000, 4000, 0, 0)
w3j = get_wigner_array(w)

function test(w_, w3j_)
    for i in 1:2000
        wigner3j_f!(w_, w3j_)
    end
end

@btime test($w, $w3j)

##

function test2(w_, w3j_)
    for i in 1:2000
        classical_wigner3j_m0!(w_, w3j_)
    end
end

@btime test2($w, $w3j)

##