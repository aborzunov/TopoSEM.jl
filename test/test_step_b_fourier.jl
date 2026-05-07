@testset "step B — Fourier integration" begin
    # Build a smooth bump h_true, derive its gradient analytically, integrate
    # and check that we recover h up to an additive constant.
    M, N = 96, 128
    dx = 1.0
    xgrid, ygrid = TopoSEM.spatial_grid(M, N, dx)
    σ = 18.0
    A = 4.0
    h_true = [A * exp(-(x^2 + y^2) / (2σ^2)) for y in ygrid, x in xgrid]
    γx = [-(x / σ^2) * h_true[i, j] for (i, y) in enumerate(ygrid), (j, x) in enumerate(xgrid)]
    γy = [-(y / σ^2) * h_true[i, j] for (i, y) in enumerate(ygrid), (j, x) in enumerate(xgrid)]
    κ  = zeros(size(h_true))   # not used when ω is small

    h_recon = TopoSEM.integrate_fourier(γx, γy, κ; dx = dx, dy = dx, ω = dx,
                                         padding = :mirror)
    # Compare after subtracting means (gauge-free DC mode)
    diff = (h_recon .- (sum(h_recon) / length(h_recon))) .-
           (h_true  .- (sum(h_true)  / length(h_true)))
    rms = sqrt(sum(abs2, diff) / length(diff))
    @test rms / maximum(abs, h_true) < 0.02

    # :none padding has visible periodic-wrap artefacts but should still be sane
    h_none = TopoSEM.integrate_fourier(γx, γy, κ; dx = dx, dy = dx, ω = dx,
                                        padding = :none)
    @test all(isfinite, h_none)
end
