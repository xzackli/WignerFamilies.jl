
using WignerFamilies
using HalfIntegers
# using OffsetArrays

# wsv = WignerSymbolVector(zeros(length(HalfInt(3/2):HalfInt(7/2))), HalfInt(3/2):HalfInt(7/2))


w = WignerF( HalfInt(5/2), 5, HalfInt(1/2), HalfInt(-1) )
m₁ = -w.m₂ - w.m₃

w3j = get_wigner_array(w)
wigner3j_f!(w, w3j)

using WignerSymbols
# @test wigner3j(Float64, HalfInt(5/2), w.j₂, w.j₃, m₁, w.m₂ ) == 

##