@testset "end-to-end on synthetic Gaussian" begin
    Random.seed!(0)
    M, N = 96, 128
    dx = 1.0
    xgrid, ygrid = TopoSEM.spatial_grid(M, N, dx)
    σ = 20.0
    A = 6.0
    h_true = [A * exp(-(x^2 + y^2) / (2σ^2)) for y in ygrid, x in xgrid]
    gx, gy = TopoSEM.gradient_centered(h_true, dx)

    angles = [0.0, π/2, π, 3π/2]
    images = Matrix{Float64}[]
    σ_noise = 0.01
    for θ in angles
        cθ, sθ = cos(θ), sin(θ)
        img = [128.0 + 1.0 * (cθ*gx[i,j] + sθ*gy[i,j]) + σ_noise * randn()
                for i in 1:M, j in 1:N]
        push!(images, img)
    end

    result = topo_sem(images, angles; pixel_size = dx, max_iter = 25,
                       tol = 1e-5, verbose = false)

    h = result.h
    hc = h .- (sum(h) / length(h))
    htc = h_true .- (sum(h_true) / length(h_true))
    rel_rms = sqrt(sum(abs2, hc .- htc) / length(h)) / maximum(abs, htc)
    @test rel_rms < 0.05
    @test result.iterations ≥ 2
    # detector vectors should align with their angle
    for d in result.detectors
        @test abs(d.p₂[2]) < 0.05            # transverse component small
        @test 0.5 < d.p₂[1] < 1.5            # n-component close to gauge
    end
end
