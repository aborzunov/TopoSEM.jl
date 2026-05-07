using Test
using TopoSEM
using LinearAlgebra
using Random

@testset "TopoSEM.jl" begin
    include("test_derivatives.jl")
    include("test_step_b_fourier.jl")
    include("test_step_c_calibration.jl")
    include("test_synthetic_gaussian.jl")
    include("test_csv_io.jl")
end
