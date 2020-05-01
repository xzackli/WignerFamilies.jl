var documenterSearchIndex = {"docs":
[{"location":"#","page":"Home","title":"Home","text":"CurrentModule = WignerFamilies","category":"page"},{"location":"#WignerFamilies-1","page":"Home","title":"WignerFamilies","text":"","category":"section"},{"location":"#","page":"Home","title":"Home","text":"","category":"page"},{"location":"#","page":"Home","title":"Home","text":"Modules = [WignerFamilies]","category":"page"},{"location":"#WignerFamilies.get_wigner_array-Union{Tuple{AbstractWigner{T}}, Tuple{T}} where T","page":"Home","title":"WignerFamilies.get_wigner_array","text":"get_wigner_array(w::AbstractWigner{T}) where {T}\n\nUtility function for getting an OffsetArray with indices from jₘᵢₙ to jₘₐₓ.\n\nArguments\n\nw::AbstractWigner{T}: contains the quantum numbers and dispatches on the kind of symbol\n\nReturns\n\nOffsetArray{T}: an array for wigner symbols\n\n\n\n\n\n","category":"method"},{"location":"#WignerFamilies.wigner3j_f-Union{Tuple{T}, Tuple{Type{T},Any,Any,Any,Any}} where T<:Real","page":"Home","title":"WignerFamilies.wigner3j_f","text":"wigner3j_f(::Type{T}, j₂, j₃, m₂, m₃) where {T<:Real}\n\nComputes all allowed j₁ given fixed j₂, j₃, m₂, m₃, m₁=-m₂-m₃. \n\nArguments\n\nT::Type{<:Real}: output array type\nj₂: quantum number\nj₃: quantum number\nm₂: quantum number\nm₃: quantum number\n\nReturns\n\nTuple{Vector{Int}, Vector{T}}: j₁ values and wigner symbols\n\n\n\n\n\n","category":"method"},{"location":"#WignerFamilies.classical_wigner3j_m0-Union{Tuple{T}, Tuple{Type{T},Any,Any,Any,Any}} where T<:Real","page":"Home","title":"WignerFamilies.classical_wigner3j_m0","text":"classical_wigner3j_m0(::Type{T}, j₂, j₃, m₂, m₃) where {T<:Real}\n\nComputes all allowed j₁ given fixed j₂, j₃, m₁, m₂, m₃, subject to m₁ + m₂ + m₃ = 0. This  applies the classical three-term recurrence relation and iterates two at a time, since in  this special case all symbols with odd ∑jᵢ are zero. Unlike other Wigner symbols, this  special case requires iterating outwards, as one must recur towards increasing |fⱼ| for  stability.\n\nArguments\n\nT::Type{<:Real}: output array type\nj₂::Tn: quantum number\nj₃::Tn: quantum number\nm₂::Tn: quantum number\n\nReturns\n\nTuple{Vector{Int}, Vector{T}}: j₁ values and wigner symbols\n\n\n\n\n\n","category":"method"},{"location":"#WignerFamilies.f_to_min_m0!-Union{Tuple{T}, Tuple{AbstractWignerF{T},Int64,AbstractArray{T,1}}} where T","page":"Home","title":"WignerFamilies.f_to_min_m0!","text":"Special case iteration for mᵢ=0.\n\n\n\n\n\n","category":"method"},{"location":"#WignerFamilies.rψ!-Union{Tuple{T}, Tuple{AbstractWigner{T},Integer,AbstractArray{T,1}}} where T","page":"Home","title":"WignerFamilies.rψ!","text":"rψ!(w::AbstractWigner{T}, n::Integer, iterates::AbstractVector{T}) where T\n\nBackward recurrence scheme defined by equation 2 in Luscombe and Luban 1998. Iteratively generates the ratio rψ(n) = ψ(n) / ψ(n+1) in the iterates vector.\n\nArguments\n\nw::AbstractWigner{T}: contains the quantum numbers and dispatches on the kind of symbol\nnmid::Integer: current index of the recurrence\nψ::AbstractVector{T}: store the values of rψ here during recursion.\n\nReturns\n\nstop::Int: the index the iteration stopped\n\n\n\n\n\n","category":"method"},{"location":"#WignerFamilies.swap_triangular-Tuple{Any}","page":"Home","title":"WignerFamilies.swap_triangular","text":"Evens out an array which scales linearly with difficulty by swapping elements such that [1,2,3,4,5,6] is mapped to [1,6,2,5,3,4].\n\n\n\n\n\n","category":"method"},{"location":"#WignerFamilies.ψauxplus!-Union{Tuple{T}, Tuple{AbstractWignerF{T},Int64,Int64,AbstractArray{T,1}}} where T","page":"Home","title":"WignerFamilies.ψauxplus!","text":"Three-term recurrence relation for the classical region.\n\n\n\n\n\n","category":"method"}]
}
