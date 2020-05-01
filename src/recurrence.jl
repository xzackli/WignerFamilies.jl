
abstract type AbstractWigner{T<:Real} end
abstract type AbstractWignerF{T} <: AbstractWigner{T} end
abstract type AbstractWignerH{T} <: AbstractWigner{T} end

struct WignerF{T} <: AbstractWignerF{T}
    j₂::T
    j₃::T
    m₂::T
    m₃::T
    nₘᵢₙ::Int
    nₘₐₓ::Int
end
function WignerF(::Type{T}, j₂, j₃, m₂, m₃) where {T<:Real}
    WignerF{T}(j₂, j₃, m₂, m₃, max(abs(j₂ - j₃), abs(m₂ + m₃)), j₂ + j₃)
end
WignerF(j₂, j₃, m₂, m₃) = WignerF(Float64, j₂, j₃, m₂, m₃)


"""
    get_wigner_array(w::AbstractWigner{T}) where {T}

Utility function for getting an OffsetArray with indices from jₘᵢₙ to jₘₐₓ.

# Arguments
- `w::AbstractWigner{T}`: contains the quantum numbers and dispatches on the kind of symbol

# Returns
- `OffsetArray{T}`: an array for wigner symbols
"""
get_wigner_array(w::AbstractWigner{T}) where {T} = OffsetArray(
    zeros(T, length(w.nₘᵢₙ:w.nₘₐₓ)), w.nₘᵢₙ:w.nₘₐₓ)


"""
    rψ!(w::AbstractWigner{T}, n::Integer, iterates::AbstractVector{T}) where T

Backward recurrence scheme defined by equation 2 in Luscombe and Luban 1998. Iteratively
generates the ratio rψ(n) = ψ(n) / ψ(n+1) in the `iterates` vector.

# Arguments
- `w::AbstractWigner{T}`: contains the quantum numbers and dispatches on the kind of symbol
- `nmid::Integer`: current index of the recurrence
- `ψ::AbstractVector{T}`: store the values of rψ here during recursion.

# Returns
- `stop::Int`: the index the iteration stopped
"""
function rψ!(w::AbstractWigner{T}, nmid::Integer, ψ::AbstractVector{T}) where T
    for n in w.nₘₐₓ:-1:nmid
        if n == w.nₘₐₓ
            ψ[n] = -Zψ(w, n) / Yψ(w, n)
        else
            ψ[n] = -Zψ(w, n) / (Yψ(w, n) + Xψ(w, n) * ψ[n+1])
        end
        if !isfinite(ψ[n])
            ψ[n] = zero(T)
            return n
        end
        if abs(ψ[n]) ≥ 1.0 
            return n
        end
    end
    return nmid
end

function sψ!(w::AbstractWigner{T}, nmid::Integer, ψ::AbstractVector{T}) where T
    for n in w.nₘᵢₙ:nmid
        if n == w.nₘᵢₙ       
            ψ[n] = -Xψ(w, n) / Yψ(w, n)
        else
            ψ[n] = -Xψ(w, n) / (Yψ(w, n) + Zψ(w, n) * ψ[n-1])
        end
        if !isfinite(ψ[n])
            ψ[n] = zero(T)
            return n
        end
        if abs(ψ[n]) ≥ one(T)
            return n
        end
    end
    return nmid
end


"""
Three-term recurrence relation for the classical region.
"""
function ψauxplus!(w::AbstractWignerF{T}, 
                   n₋::Int, nc::Int, ψ::AbstractVector{T}) where {T}
    start_index = n₋
    if n₋ == w.nₘᵢₙ
        ψ[w.nₘᵢₙ] = one(T)  # set an arbitrary value for this
        if iszero(w.nₘᵢₙ)
            ψ[w.nₘᵢₙ+1] = -(w.m₃-w.m₂+2B(w,w.nₘᵢₙ)) / A(w,1)
        else
            ψ[w.nₘᵢₙ+1] = -(B(w,w.nₘᵢₙ)) / (w.nₘᵢₙ* A(w, w.nₘᵢₙ+1))
        end
        start_index = w.nₘᵢₙ+1
    end
    
    for n in start_index:(nc-1)
        Xn = Xψ(w, n)
        if iszero(Xn)
            ψ[n+1] = zero(T)
        else
            ψ[n+1] = -(Yψ(w, n) * ψ[n] + Zψ(w, n) * ψ[n - 1]) / Xn
        end
    end
end



A(w::AbstractWignerF, j) = sqrt((j^2 - (w.j₂ - w.j₃)^2) * 
    ((w.j₂ + w.j₃ + 1)^2 - j^2) * (j^2 - (w.m₂ + w.m₃)^2))
B(w::AbstractWignerF, j) = (2 * j + 1) * ((w.m₂ + w.m₃) * (w.j₂ * (w.j₂ + 1) - 
    w.j₃ * (w.j₃ + 1)) - (w.m₂ - w.m₃) * j * (j + 1))
Xψ(w::AbstractWignerF, j) = j * A(w, j+1)
Yψ(w::AbstractWignerF, j) = B(w, j)
Zψ(w::AbstractWignerF, j) = (j + 1) * A(w, j)
function normalization(w::AbstractWignerF{T}, ψ₀::AbstractVector{T}) where T
    norm = zero(T)
    for j in eachindex(ψ₀)
        norm += (2j + 1) * ψ₀[j]^2
    end
    return sqrt(abs(norm))
end
f_jmax_sgn(w::AbstractWignerF) = iseven(Int(w.j₂ - w.j₃ + w.m₂ + w.m₃)) ? 1 : -1

struct WignerH{T} <: AbstractWignerH{T}
    j₂::T
    j₃::T
    l₁::T
    l₂::T
    l₃::T
    nₘᵢₙ::T
    nₘₐₓ::T

    WignerH(::Type{T}, j₂, j₃, l₁, l₂, l₃) where {T} = new{T}(
        j₂, j₃, l₁, l₂, l₃, max(abs(j₂ - j₃), abs(l₂ - l₃)), min(j₂ + j₃, l₂ + l₃))
end

E(w::AbstractWignerH, j) = sqrt((j^2 - (w.j₂ - w.j₃)^2) * ((w.j₂ + w.j₃ + 1)^2 - j^2) * 
    (j^2 - (w.l₂ - w.l₃)^2) * ((w.l₂ + w.l₃ + 1)^2 - j^2))
F(w::AbstractWignerH, j) = (2j+1) * (
    j * (j + 1)  * (-j * (j+1) + w.j₂ * (w.j₂ + 1) + w.j₃ * (w.j₃ + 1) - 2 * w.l₁ * (w.l₁ + 1)) + 
    w.l₂ * (w.l₂ + 1) * (j * (j+1) + w.j₂ * (w.j₂ + 1) - w.j₃ * (w.j₃ + 1)) +
    w.l₃ * (w.l₃ + 1) * (j * (j+1) - w.j₂ * (w.j₂ + 1) + w.j₃ * (w.j₃ + 1))
)
Xψ(w::AbstractWignerH, j) = j * E(w, j+1)
Yψ(w::AbstractWignerH, j) = F(w, j)
Zψ(w::AbstractWignerH, j)= (j+1) * E(w, j)

"""
    nonclassical_wigner3j(::Type{T}, j₂, j₃, m₂, m₃) where {T<:Real}

Computes all allowed j₁ given fixed j₂, j₃, m₂, m₃, m₁=-m₂-m₃. This only is guarantted to
work in non-classical regions.

# Arguments
- `T::Type{<:Real}`: output array type
- `j₂`: quantum number
- `j₃`: quantum number
- `m₂`: quantum number
- `m₃`: quantum number

# Returns
- `Tuple{Vector{Int}, Vector{T}}`: j₁ values and wigner symbols
"""
function nonclassical_wigner3j(::Type{T}, j₂, j₃, m₂, m₃) where {T<:Real}
    w = WignerF(T, j₂, j₃, m₂, m₃)
    w3j = get_wigner_array(w)
    nonclassical_wigner3j!(w, w3j)
    return collect(w.nₘᵢₙ:w.nₘₐₓ), w3j
end
nonclassical_wigner3j(j₂, j₃, m₂, m₃) = nonclassical_wigner3j(Float64, j₂, j₃, m₂, m₃)
function nonclassical_wigner3j!(w::AbstractWignerF{T}, w3j::AbstractVector{T}) where T
    nmid = Int(ceil((w.nₘᵢₙ + w.nₘₐₓ)/2))
    rψ!(w, nmid, w3j)
    sψ!(w, nmid, w3j)
    for i in (nmid+1):w.nₘₐₓ  # product the ratios upwards
        w3j[i] *= w3j[i-1]
    end
    for i in (nmid-1):-1:w.nₘᵢₙ  # product the ratios downwards
        w3j[i] *= w3j[i+1]
    end
    norm = normalization(w, w3j)
    w3j ./= norm
    if sign(w3j[w.nₘₐₓ]) != f_jmax_sgn(w)
        w3j .*= -1
    end
end


"""
Special case iteration for mᵢ=0.
"""
function f_to_min_m0!(w::AbstractWignerF{T}, 
                      nmid::Int, ψ::AbstractVector{T}) where {T}
    @assert iseven(nmid)
    ψ[nmid] = one(T)
    ψ[nmid+1] = zero(T)
    for n in (nmid+1):2:(w.nₘₐₓ-1)
        Xn = Xψ(w, n)
        if Xn == 0
            ψ[n+1] = 0
        else
            ψ[n+1] = -(Yψ(w, n) * ψ[n] + Zψ(w, n) * ψ[n - 1]) / Xn
        end
    end
end

function f_to_max_m0!(w::AbstractWignerF{T}, nmid::Int, ψ::AbstractVector{T}) where {T}
    @assert iseven(nmid)
    ψ[nmid] = one(T)
    ψ[nmid-1] = zero(T)
    for n in (nmid-1):-2:(w.nₘᵢₙ+1)
        Zn = Zψ(w, n)
        if Zn == 0
            ψ[n-1] = 0
        else
            ψ[n-1] = -(Yψ(w, n) * ψ[n] + Xψ(w, n) * ψ[n + 1]) / Zn
        end
    end
end


"""
    classical_wigner3j_m0(::Type{T}, j₂, j₃, m₂, m₃) where {T<:Real}

Computes all allowed j₁ given fixed j₂, j₃, m₁, m₂, m₃, subject to m₁ + m₂ + m₃ = 0. This 
applies the classical three-term recurrence relation and iterates two at a time, since in 
this special case all symbols with odd ∑jᵢ are zero. Unlike other Wigner symbols, this 
special case requires iterating outwards, as one must recur towards increasing |fⱼ| for 
stability.

# Arguments
- `T::Type{<:Real}`: output array type
- `j₂::Tn`: quantum number
- `j₃::Tn`: quantum number
- `m₂::Tn`: quantum number

# Returns
- `Tuple{Vector{Int}, Vector{T}}`: j₁ values and wigner symbols
"""
function classical_wigner3j_m0(::Type{T}, j₂, j₃, m₂, m₃) where {T<:Real}
    w = WignerF(j₂, j₃, m₂, m₃)
    w3j = get_wigner_array(w)
    classical_wigner3j_m0!(w, w3j)
    return collect(w.nₘᵢₙ:w.nₘₐₓ), w3j
end

function classical_wigner3j_m0!(w::AbstractWignerF{T}, w3j::AbstractVector{T}) where T
    nmid = Int( (w.nₘᵢₙ + w.nₘₐₓ) / 2 )
    nmid = iseven(nmid) ? nmid : nmid + 1  # ensure start index is even
    f_to_min_m0!(w, nmid, w3j)
    f_to_max_m0!(w, nmid, w3j)
    norm = normalization(w, w3j)
    w3j ./= norm
    if sign(w3j[w.nₘₐₓ]) != f_jmax_sgn(w)
        w3j .*= -1
    end
end

