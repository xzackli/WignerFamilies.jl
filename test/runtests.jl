using Test
using WignerFamilies
using HalfIntegers
import WignerSymbols  # only a test-dep, for comparisons

@testset "f: nonclassical" begin
    j₂ = 100
    j₃ = 60
    m₂ = 70
    m₃ = -55
    m₁ = -m₂ - m₃
    w3j = wigner3j_f(Float64, j₂, j₃, m₂, m₃)
    j_array = collect(eachindex(w3j))
    for j in j_array
        @test w3j[j] ≈ WignerSymbols.wigner3j(Float64, j, j₂, j₃, m₁, m₂)
    end
end

@testset "f: ∑mᵢ = 0 special case" begin
    j₂, j₃, m₂, m₃ = Int128.((5000, 5002, 0, 0))
    w3j = wigner3j_f(Float64, j₂, j₃, m₂, m₃)
    j_array = collect(eachindex(w3j))
    for j in [j_array[i] for i in [1, 2, 4, 6, lastindex(j_array)]]
        @test w3j[j] ≈ Float64(WignerSymbols.wigner3j(BigFloat, j, j₂, j₃, -m₂-m₃, m₂))
    end
end

## bottom of page 14 of Raynal et al. On the zeros of 3j coefficients: polynomial degree 
# versus recurrence order is an example of sequence of nontrivial zeros
@testset "f: nontrivial zeros" begin
    tol = eps()
    for (X₀, Y₀) in ((43, 25), (109, 63))
        j₂, j₃, m₂, m₃ = Int128.((5000, 5002, 0, 0))
        X(n) = iszero(n) ? X₀ : 7X(n-1) + 12Y(n-1)
        Y(n) = iszero(n) ? Y₀ : 4(n-1) + 7Y(n-1)

        for n_seq in 0:1
            j₂, j₃, m₂, m₃ = (2X(n_seq)+1)/6, (X(n_seq)+2)/6, HalfInt(3/2), HalfInt(3/2)
            w3j = wigner3j_f(Float64, j₂, j₃, m₂, m₃)
            j_array = collect(eachindex(w3j))
            for j in j_array
                @test (abs(w3j[j] - 
                    Float64(WignerSymbols.wigner3j(BigFloat, j, j₂, j₃, -m₂-m₃, m₂))) < tol)
            end
        end
    end
end

## tests for the Rasch and Yu c index
function confirm_symmetry(maxj)
    j₁, j₂ = rand(0:maxj, 2)
    j₃ = rand(abs(j₁ - j₂):(j₁ + j₂))
    if isodd(j₁ + j₂ + j₃)
        j₃ += 1
    end
    m₁, m₂, m₃ = 0, 0, 0
    c1 = rasch_yu_index(Int128, j₁, j₂, j₃, m₁, m₂, m₃)
    c2 = rasch_yu_index(Int128, 
        j₁, (j₂ + j₃ - m₁)/2, (j₂ + j₃ + m₁)/2, 
        j₃ - j₂, (j₂ - j₃ - m₁)/2 - m₃, (j₂ - j₃ + m₁)/2 + m₃)

    c1, c2
end

@testset "Rasch & Yu c: Regge symmetry" begin
    for i in 1:1000
        c1, c2 = confirm_symmetry(1000)
        @test c1 == c2
    end
end
