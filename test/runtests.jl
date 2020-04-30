using WignerFamilies
import WignerSymbols
using Test

@testset "f: nonclassical" begin
    j₂ = 100
    j₃ = 60
    m₂ = 70
    m₃ = -55
    m₁ = -m₂ - m₃
    j_array, w3j = nonclassical_wigner3j(Float64, j₂, j₃, m₂, m₃)
    for j in j_array
        @test w3j[j] ≈ WignerSymbols.wigner3j(Float64, j, j₂, j₃, m₁, m₂)
    end
end

@testset "f: ∑mᵢ = 0" begin
    j₂, j₃, m₂, m₃ = Int128.((5000, 5002, 0, 0))
    j_array, w3j = WignerFamilies.classical_wigner3j_m0(Float64, j₂, j₃, m₂, m₃)
    for j in [j_array[1], j_array[2], j_array[3], j_array[end]]
        @test w3j[j] ≈ Float64(WignerSymbols.wigner3j(BigFloat, j, j₂, j₃, -m₂-m₃, m₂))
    end
end

