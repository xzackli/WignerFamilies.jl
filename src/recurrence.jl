using HalfIntegers
using OffsetArrays

abstract type AbstractWigner{T<:Union{Half{Integer}, Integer}} end
abstract type AbstractWignerF{T} <: AbstractWigner{T} end
abstract type AbstractWignerH{T} <: AbstractWigner{T} end

struct WignerF{T} <: AbstractWignerF{T}
    j₂::T
    j₃::T
    m₂::T
    m₃::T
    nₘᵢₙ::T
    nₘₐₓ::T
end
WignerF(j₂::T, j₃::T, m₂::T, m₃::T) where T = WignerF{T}(
    j₂, j₃, m₂, m₃, max(abs(j₂ - j₃), abs(m₂ + m₃)), j₂ + j₃)

function rψ!(w::AbstractWigner, n::Integer, iterates::AbstractVector{<:Real})
    if n == w.nₘₐₓ
        iterates[n] = -Zψ(w, n) / Yψ(w, n)
        return iterates[n]
    end
    @assert n < w.nₘₐₓ
    t = rψ!(w, n + 1, iterates)
    iterates[n] = -Zψ(w, n) / (Yψ(w, n) + Xψ(w, n) * t)
    return iterates[n] 
end

function sψ!(w::AbstractWigner, n::Integer, iterates::AbstractVector{<:Real})
    if n == w.nₘᵢₙ
        iterates[n] = -Xψ(w, n) / Yψ(w, n)
        return iterates[n]
    end
    @assert n > w.nₘᵢₙ
    t = sψ!(w, n - 1, iterates)
    iterates[n] = -Xψ(w, n) / (Yψ(w, n) + Zψ(w, n) * t)
    return iterates[n]
end

A(w::AbstractWignerF, j) = sqrt((j^2 - (w.j₂ - w.j₃)^2) * 
    ((w.j₂ + w.j₃ + 1)^2 - j^2) * (j^2 - (w.m₂ + w.m₃)^2))
B(w::AbstractWignerF, j) = (2 * j + 1) * ((w.m₂ + w.m₃) * (w.j₂ * (w.j₂ + 1) - 
    w.j₃ * (w.j₃ + 1)) - (w.m₂ - w.m₃) * j * (j + 1))
Xψ(w::AbstractWignerF, j) = j * A(w, j+1)
Yψ(w::AbstractWignerF, j) = B(w, j)
Zψ(w::AbstractWignerF, j) = (j + 1) * A(w, j)
function normalization(w::AbstractWignerF, ψ₀::AbstractVector{<:Real})
    norm = big(0.0)
    for j in eachindex(ψ₀)
        norm += (2j + 1) * ψ₀[j]^2
    end
    return sqrt(abs(norm))
end


struct WignerH{T} <: AbstractWignerH{T}
    j₂::T
    j₃::T
    l₁::T
    l₂::T
    l₃::T
    nₘᵢₙ::T
    nₘₐₓ::T

    WignerH(j₂::T, j₃::T, l₁::T, l₂::T, l₃::T) where T = new{T}(
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
    nonclassical_wigner3j(T::Type{<:Real}, j₂::Tn, j₃::Tn, m₂::Tn, m₃::Tn) where {T, Tn}

Computes all allowed j₁ given fixed j₂, j₃, m₂, m₃, m₁=-m₂-m₃. This only is guarantted to
work in non-classical regions.

# Arguments
- `T::Type{<:Real}`: output array type
- `j₂::Tn`: quantum number
- `j₃::Tn`: quantum number
- `m₂::Tn`: quantum number
- `m₃::Tn`: quantum number

# Returns
- `Tuple{Vector{Int}, Vector{T}}`: j₁ values and wigner symbols
"""
function nonclassical_wigner3j(T::Type{<:Real}, 
                               j₂::Tn, j₃::Tn, m₂::Tn, m₃::Tn) where {Tn}
    w = WignerF(j₂, j₃, m₂, m₃)
    nmid = Int(ceil((w.nₘᵢₙ + w.nₘₐₓ)/2))
    w3j = OffsetArray(zeros(T, length(w.nₘᵢₙ:w.nₘₐₓ)), w.nₘᵢₙ:w.nₘₐₓ)
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
    fjmax_sgn = iseven(w.j₂ + w.j₃ + w.m₂ + w.m₃) ? 1 : -1
    if sign(w3j[w.nₘₐₓ]) != fjmax_sgn
        w3j .*= -1
    end
    return collect(w.nₘᵢₙ:w.nₘₐₓ), w3j
end


function f_to_min_m0!(w::AbstractWignerF, 
                      nmid::Int, ψ::AbstractVector{T}) where {T<:Real}
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

function f_to_max_m0!(w::AbstractWignerF, nmid::Int, 
                      ψ::AbstractVector{T}) where {T<:Real}
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
    classical_wigner3j_m0(T::Type{<:Real}, j₂::Tn, j₃::Tn, m₂::Tn, m₃::Tn) where {T, Tn}

Computes all allowed j₁ given fixed j₂, j₃, m₁ + m₂ + m₃ = 0. This applies the classical
three-term recurrence relation and iterates two at a time, since all odd ∑jᵢ are zero.
Unlike other Wigner symbols, this special case requires iterating outwards, as one must
recur towards increasing |fⱼ| for stability.

# Arguments
- `T::Type{<:Real}`: output array type
- `j₂::Tn`: quantum number
- `j₃::Tn`: quantum number
- `m₂::Tn`: quantum number

# Returns
- `Tuple{Vector{Int}, Vector{T}}`: j₁ values and wigner symbols
"""
function classical_wigner3j_m0(T::Type{<:Real}, 
                               j₂::Tn, j₃::Tn, m₂::Tn, m₃::Tn) where {Tn}
    
    w = WignerF(j₂, j₃, m₂, m₃)
    nmid = Int( (w.nₘᵢₙ + w.nₘₐₓ) / 2 )
    nmid = iseven(nmid) ? nmid : nmid + 1
    j_array = collect(Int(w.nₘᵢₙ):Int(w.nₘₐₓ))
    w3j = OffsetArray(zeros(T, length(j_array)), j_array[1]:j_array[end])
    f_to_min_m0!(w, nmid, w3j)
    f_to_max_m0!(w, nmid, w3j)
    norm = normalization(w, w3j)
    w3j ./= norm
    fjmax_sgn = iseven(w.j₂ + w.j₃ + w.m₂ + w.m₃) ? 1 : -1
    if sign(w3j[w.nₘₐₓ]) != fjmax_sgn
        w3j .*= -1
    end
    return j_array, w3j
end

