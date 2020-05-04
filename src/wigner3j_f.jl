
# the result arrays are indexed by type Ti and have type T
struct WignerF{T<:Real, Ti} <: AbstractWignerF{T, Ti}
    j₂::T
    j₃::T
    m₂::T
    m₃::T
    nₘᵢₙ::Ti
    nₘₐₓ::Ti
end
function WignerF(::Type{T}, j₂::Ti, j₃::Ti, m₂::Ti, m₃::Ti) where {T<:Real, Ti}
    WignerF{T, Ti}(j₂, j₃, m₂, m₃, max(abs(j₂ - j₃), abs(m₂ + m₃)), j₂ + j₃)
end
function WignerF(::Type{T}, j₂, j₃, m₂, m₃) where {T}
    nₘᵢₙ = max(abs(j₂ - j₃), abs(m₂ + m₃))
    nₘₐₓ = j₂ + j₃
    # if isinteger(nₘᵢₙ) && isinteger(nₘₐₓ)
    #     return WignerF{T, Int}(j₂, j₃, m₂, m₃, nₘᵢₙ, nₘₐₓ)
    # else
    return WignerF{T, HalfInt}(j₂, j₃, m₂, m₃, nₘᵢₙ, nₘₐₓ)  # better to just be type stable
    # end
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
get_wigner_array(w::AbstractWigner{T,Ti}) where {T,Ti} = (
    WignerSymbolVector(zeros(T, length(w.nₘᵢₙ:w.nₘₐₓ)), w.nₘᵢₙ:w.nₘₐₓ))


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
function rψ!(w::AbstractWigner{T,Ti}, nmid::Ti, ψ::AbstractVector{T}) where {T,Ti}
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

function sψ!(w::AbstractWigner{T,Ti}, nmid::Ti, ψ::AbstractVector{T}) where {T,Ti}
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
function ψauxplus!(w::AbstractWignerF{T, Ti}, 
                   n₋::Ti, nc::Ti, ψ::AbstractVector{T}) where {T, Ti}
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


function selection_rules(w::AbstractWignerF{T}) where T
    small = 10 * eps(T)
    return (
        (abs(w.m₂) < w.j₂ + small) &&
        (abs(w.m₃) < w.j₃ + small) &&
        isinteger(w.m₂ - w.j₂) &&
        isinteger(w.m₃ - w.j₃)
    )
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

"""
    wigner3j_f(::Type{T}, j₂, j₃, m₂, m₃) where {T<:Real}

Computes all allowed j₁ given fixed j₂, j₃, m₂, m₃, m₁=-m₂-m₃. 

# Arguments
- `T::Type{<:Real}`: output array type
- `j₂`: quantum number
- `j₃`: quantum number
- `m₂`: quantum number
- `m₃`: quantum number

# Returns
- `Tuple{Vector{Int}, Vector{T}}`: j₁ values and wigner symbols
"""
function wigner3j_f(::Type{T}, j₂, j₃, m₂, m₃) where {T<:Real}
    w = WignerF(T, j₂, j₃, m₂, m₃)  # establish the problem variables and method
    if !selection_rules(w)
        # throw(DomainError("Bad combination of inputs."))
        return OffsetArray(T[], 1:0)
    end
    w3j = get_wigner_array(w)  # construct an OffsetArray with the right indices and type
    wigner3j_f!(w, w3j)  # call the mutating function
    return w3j
end
wigner3j_f(j₂, j₃, m₂, m₃) = wigner3j_f(Float64, j₂, j₃, m₂, m₃)

function wigner3j_f!(w::AbstractWignerF{T,Ti}, w3j::AbstractVector{T}) where {T,Ti}

    # special case, if it's just one nontrivial result
    if length(w3j) == 1
        w3j[w.nₘᵢₙ] = f_jmax_sgn(w) / sqrt(w.nₘᵢₙ + w.j₂ + w.j₃+1)
        return
    end

    nmid = Int( (w.nₘᵢₙ + w.nₘₐₓ) / 2 )
    fill!(w3j.symbols, zero(T))

    # special case that performs an outwards classical solution if m₁ = m₂ = m₃ = 0
    # and skips the symbols with odd ∑jᵢ since those are zero.
    if iszero(w.m₂) && iszero(w.m₃)
        nmid_even = iseven(nmid) ? nmid : nmid + 1  # ensure start index is even
        if w.nₘᵢₙ < nmid_even < w.nₘₐₓ  # check to make sure it's not the ends
            return classical_wigner3j_m0!(w, w3j)
        end
    end

    # otherwise attempt the general solution
    nmid = Ti(ceil((w.nₘᵢₙ + w.nₘₐₓ)/2))
    # attempt the non-classical two term nonlinear recurrences
    n₊ = rψ!(w, nmid, w3j) 
    n₋ = sψ!(w, nmid, w3j)
    for k in (n₊+1):w.nₘₐₓ  # product the ratios upwards
        w3j[k] *= w3j[k-1]
    end
    for k in (n₋-1):-1:w.nₘᵢₙ  # product the ratios downwards
        w3j[k] *= w3j[k+1]
    end
    # now use the three-term classical recurrence in the middle
    ψn₋ = w3j[n₋]
    ψn₊ = w3j[n₊]
    ψauxplus!(w, n₋, n₊, w3j)
    # rescale the three parts
    if n₋ > w.nₘᵢₙ
        scale_middle = ψn₋ / w3j[n₋]
        for j in w.nₘᵢₙ:(n₋+1)
            w3j[j] *= scale_middle
        end
    end
    if n₊ < w.nₘₐₓ
        scale_end = w3j[n₊] / ψn₊
        for j in (n₊+1):w.nₘₐₓ
            w3j[j] *= scale_end
        end
    end
    # normalize the results
    norm = normalization(w, w3j)
    w3j.symbols ./= norm
    if sign(w3j[w.nₘₐₓ]) != f_jmax_sgn(w)
        w3j.symbols .*= -1
    end
end


"""
Special case iteration for mᵢ=0. In this case, j₁ must be an integer too, so the index
type is restricted to Int.
"""
function f_to_min_m0!(w::AbstractWignerF{T, Int}, 
                      nmid::Int, ψ::AbstractVector{T}) where {T}
    @assert iseven(nmid)
    ψ[nmid] = one(T)
    ψ[nmid+1] = zero(T)
    for n in (nmid+1):2:(w.nₘₐₓ-1)
        Xn = Xψ(w, n)
        @assert abs(Xn) > zero(T)
        ψ[n+1] = -(Yψ(w, n) * ψ[n] + Zψ(w, n) * ψ[n - 1]) / Xn
    end
end

function f_to_max_m0!(w::AbstractWignerF{T, Int}, 
                      nmid::Int, ψ::AbstractVector{T}) where {T}
    @assert iseven(nmid)
    ψ[nmid] = one(T)
    ψ[nmid-1] = zero(T)
    for n in (nmid-1):-2:(w.nₘᵢₙ+1)
        Zn = Zψ(w, n)
        @assert abs(Zn) > zero(T)
        ψ[n-1] = -(Yψ(w, n) * ψ[n] + Xψ(w, n) * ψ[n + 1]) / Zn
    end
end


"""
    classical_wigner3j_m0!(::Type{T}, j₂, j₃, m₂, m₃) where {T<:Real}

Computes all allowed j₁ given fixed j₂, j₃, m₁, m₂, m₃, subject to m₁ = m₂ = m₃ = 0. This 
applies the classical three-term recurrence relation and iterates two at a time, since in 
this special case all symbols with odd ∑jᵢ are zero. Unlike other Wigner symbols, this 
special case requires iterating outwards, as one must recur towards increasing |fⱼ| for 
stability.

# Arguments
- `T::Type{<:Real}`: output array type
- `j₂::Tn`: quantum number
- `j₃::Tn`: quantum number
- `m₂::Tn`: quantum number
"""
function classical_wigner3j_m0!(w::AbstractWignerF{T,Int}, w3j::AbstractVector{T}) where T
    nmid = Int( (w.nₘᵢₙ + w.nₘₐₓ) / 2 )
    nmid = iseven(nmid) ? nmid : nmid + 1  # ensure start index is even
    f_to_min_m0!(w, nmid, w3j)
    f_to_max_m0!(w, nmid, w3j)
    norm = normalization(w, w3j)
    w3j.symbols ./= norm
    if sign(w3j[w.nₘₐₓ]) != f_jmax_sgn(w)
        w3j.symbols .*= -1
    end
end


## ----- todo: 6j

# struct WignerH{T} <: AbstractWignerH{T}
#     j₂::T
#     j₃::T
#     l₁::T
#     l₂::T
#     l₃::T
#     nₘᵢₙ::T
#     nₘₐₓ::T

#     WignerH(::Type{T}, j₂, j₃, l₁, l₂, l₃) where {T} = new{T}(
#         j₂, j₃, l₁, l₂, l₃, max(abs(j₂ - j₃), abs(l₂ - l₃)), min(j₂ + j₃, l₂ + l₃))
# end

# E(w::AbstractWignerH, j) = sqrt((j^2 - (w.j₂ - w.j₃)^2) * ((w.j₂ + w.j₃ + 1)^2 - j^2) * 
#     (j^2 - (w.l₂ - w.l₃)^2) * ((w.l₂ + w.l₃ + 1)^2 - j^2))
# F(w::AbstractWignerH, j) = (2j+1) * (
#     j * (j + 1)  * (-j * (j+1) + w.j₂ * (w.j₂ + 1) + w.j₃ * (w.j₃ + 1) - 2 * w.l₁ * (w.l₁ + 1)) + 
#     w.l₂ * (w.l₂ + 1) * (j * (j+1) + w.j₂ * (w.j₂ + 1) - w.j₃ * (w.j₃ + 1)) +
#     w.l₃ * (w.l₃ + 1) * (j * (j+1) - w.j₂ * (w.j₂ + 1) + w.j₃ * (w.j₃ + 1))
# )
# Xψ(w::AbstractWignerH, j) = j * E(w, j+1)
# Yψ(w::AbstractWignerH, j) = F(w, j)
# Zψ(w::AbstractWignerH, j)= (j+1) * E(w, j)
