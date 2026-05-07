@testset "step C — calibration recovers params" begin
    # Build a known h, synthesise images with known per-detector params,
    # corrupt the params, run recalibrate! once, check recovery.
    Random.seed!(42)
    M, N = 64, 80
    dx = 1.0
    xgrid, ygrid = TopoSEM.spatial_grid(M, N, dx)
    σ = 14.0
    h = [3.0 * exp(-(x^2 + y^2) / (2σ^2)) for y in ygrid, x in xgrid]
    gx, gy = TopoSEM.gradient_centered(h, dx)
    lap    = TopoSEM.laplacian_5pt(h, dx)

    angles = [0.0, π/2, π, 3π/2]
    true_params = [(p₀ = 100.0 + 5k, p₂n = 1.0 + 0.1k, p₂t = 0.05k) for k in 1:4]

    images = Matrix{Float64}[]
    for (k, θ) in enumerate(angles)
        cθ, sθ = cos(θ), sin(θ)
        tp = true_params[k]
        img = [tp.p₀ + tp.p₂n * (cθ*gx[i,j] + sθ*gy[i,j]) +
                       tp.p₂t * (-sθ*gx[i,j] + cθ*gy[i,j])
                for i in 1:M, j in 1:N]
        push!(images, img)
    end

    # Detectors initialised with WRONG params
    detectors = [TopoSEM.DetectorModel(angles[k], images[k]; order = 1, max_val = 255.0)
                 for k in 1:4]
    weights = [ones(M, N) for _ in 1:4]

    TopoSEM.recalibrate!(detectors, images, h, gx, gy, lap, xgrid, ygrid, weights;
                         use_laplacian = false, use_wide_fov = false)

    for k in 1:4
        @test isapprox(detectors[k].p₀,    true_params[k].p₀;  atol = 1e-6)
        @test isapprox(detectors[k].p₂[1], true_params[k].p₂n; atol = 1e-6)
        @test isapprox(detectors[k].p₂[2], true_params[k].p₂t; atol = 1e-6)
    end
end
