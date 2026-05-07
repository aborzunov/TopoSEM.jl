@testset "derivatives" begin
    M, N = 96, 128
    L = 1.0
    dx = L / N
    dy = L / M
    xgrid = collect((1:N) .* dx)
    ygrid = collect((1:M) .* dy)
    h  = [sin(2π * x) * cos(2π * y) for y in ygrid, x in xgrid]
    gx_true = [ 2π * cos(2π * x) * cos(2π * y) for y in ygrid, x in xgrid]
    gy_true = [-2π * sin(2π * x) * sin(2π * y) for y in ygrid, x in xgrid]

    gx, gy = TopoSEM.gradient_centered(h, dx, dy)
    # interior error should be O(dx²)
    interior(A) = A[3:end-2, 3:end-2]
    err_x = maximum(abs, interior(gx) .- interior(gx_true))
    err_y = maximum(abs, interior(gy) .- interior(gy_true))
    @test err_x < 0.1
    @test err_y < 0.1

    # Laplacian of sin(2πx)cos(2πy) is -8π²·h
    lap = TopoSEM.laplacian_5pt(h, dx, dy)
    lap_true = -8π^2 .* h
    err_lap = maximum(abs, interior(lap) .- interior(lap_true))
    @test err_lap < 5.0   # 5-pt FD has visible error at this resolution
end
