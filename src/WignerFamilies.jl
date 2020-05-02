module WignerFamilies

using HalfIntegers
using OffsetArrays

export WignerSymbolVector
export WignerF, get_wigner_array
export wigner3j_f, wigner3j_f!
export rasch_yu_index

abstract type AbstractWigner{T, Ti} end
abstract type AbstractWignerF{T, Ti} <: AbstractWigner{T, Ti} end
# abstract type AbstractWignerH{T, Ti} <: AbstractWigner{T, Ti} end  # TODO: 6j

include("wignersymbolvector.jl")
include("rasch_yu.jl")
include("wigner3j_f.jl")

end
