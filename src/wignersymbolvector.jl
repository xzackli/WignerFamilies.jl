"""
This relatively light array wrapper allows us to access arrays with the quantum number 
as the index, which means we can just type in the equations from the papers.
"""

using Base: @propagate_inbounds
using CustomUnitRanges: filename_for_urange
include(filename_for_urange)

struct WignerSymbolVector{T, Ti, AA<:AbstractArray{T,1}} <: AbstractArray{T,1}
    symbols::AA
    offset::Ti
end

Base.parent(A::WignerSymbolVector) = A.symbols
offset(axparent::AbstractUnitRange, ax::AbstractUnitRange) = first(ax) - first(axparent)

function WignerSymbolVector(A::AbstractArray{T,1}, inds::AbstractUnitRange) where {T}
    axparent = first(axes(A))
    WignerSymbolVector(A, offset(axparent, inds))
end

Base.eachindex(::IndexLinear, A::WignerSymbolVector)   = axes(A, 1)
Base.IndexStyle(::Type{<:WignerSymbolVector}) = IndexLinear()
@inline Base.size(A::WignerSymbolVector) = size(parent(A))

Base.@propagate_inbounds function Base.getindex(
        A::WignerSymbolVector{T,Ti,AA}, hi::Ti) where {T, Ti<:HalfInteger, AA}
    return parent(A)[div(hi.twice - A.offset.twice, 2)]
end
Base.@propagate_inbounds function Base.getindex(
        A::WignerSymbolVector{T,Ti,AA}, hi::Rational) where {T, Ti<:HalfInteger, AA}
    return parent(A)[div(2 * hi - A.offset.twice, 2)]
end
@propagate_inbounds Base.getindex(
    A::WignerSymbolVector{T,Ti,AA}, i::Ti) where {T, Ti<:Int, AA}  = parent(A)[i - A.offset]


Base.@propagate_inbounds function Base.setindex!(
        A::WignerSymbolVector{T,Int,AA}, val, i::Int) where {T, AA}
    @boundscheck checkbounds(A, i)
    @inbounds parent(A)[i - A.offset] = val
    val
end
Base.@propagate_inbounds function Base.setindex!(
        A::WignerSymbolVector{T,HalfInt,AA}, val, hi::HalfInt) where {T, AA}
    @boundscheck checkbounds(A, hi)
    @inbounds parent(A)[div(hi.twice - A.offset.twice, 2)] = val
    val
end

Base.axes(A::WignerSymbolVector{T,Ti,AA}) where {T, Ti, AA} = 
    map(n->URange{Ti}(1+A.offset,n+A.offset), size(A))
Base.axes1(A::WignerSymbolVector{T,Ti,AA}) where {T, Ti, AA} = URange{Ti}(1+A.offset,size(A,1)+A.offset)

Base.dataids(A::WignerSymbolVector) = Base.dataids(parent(A))

function Base.show(io::IO, ::MIME"text/plain", a::WignerSymbolVector)
    print(io, "$(length(a))-element WignerSymbolVector (")
    println(io, eltype(a), ") over ", Base.axes1(a))
    Base.print_array(io, parent(a))
end
