"""
    gradient_centered(h, dx, dy = dx) -> (gx, gy)

Second-order centered finite differences in the interior, second-order
one-sided differences at the boundary. Convention: `x` runs along columns
(index `j`), `y` runs along rows (index `i`).

`gx[i,j] = ∂h/∂x` and `gy[i,j] = ∂h/∂y`.
"""
function gradient_centered(h::AbstractMatrix{T}, dx::Real, dy::Real = dx) where {T}
    M, N = size(h)
    gx = similar(h)
    gy = similar(h)
    inv2dx = T(1 / (2 * dx))
    inv2dy = T(1 / (2 * dy))
    @inbounds for j in 1:N, i in 1:M
        if j == 1
            gx[i,j] = (-3h[i,1] + 4h[i,2] - h[i,3]) * inv2dx
        elseif j == N
            gx[i,j] = (3h[i,N] - 4h[i,N-1] + h[i,N-2]) * inv2dx
        else
            gx[i,j] = (h[i,j+1] - h[i,j-1]) * inv2dx
        end
        if i == 1
            gy[i,j] = (-3h[1,j] + 4h[2,j] - h[3,j]) * inv2dy
        elseif i == M
            gy[i,j] = (3h[M,j] - 4h[M-1,j] + h[M-2,j]) * inv2dy
        else
            gy[i,j] = (h[i+1,j] - h[i-1,j]) * inv2dy
        end
    end
    return gx, gy
end

"""
    laplacian_5pt(h, dx, dy = dx)

Five-point Laplacian with reflection boundary conditions:
`∇²h[i,j] ≈ (h[i+1,j] + h[i-1,j])/dy² + (h[i,j+1] + h[i,j-1])/dx² - 2(1/dx²+1/dy²)·h[i,j]`.
"""
function laplacian_5pt(h::AbstractMatrix{T}, dx::Real, dy::Real = dx) where {T}
    M, N = size(h)
    out = similar(h)
    inv_dx2 = T(1 / dx^2)
    inv_dy2 = T(1 / dy^2)
    @inbounds for j in 1:N, i in 1:M
        ip = i == M ? i - 1 : i + 1   # reflect
        im = i == 1 ? i + 1 : i - 1
        jp = j == N ? j - 1 : j + 1
        jm = j == 1 ? j + 1 : j - 1
        out[i,j] = (h[ip,j] + h[im,j] - 2h[i,j]) * inv_dy2 +
                   (h[i,jp] + h[i,jm] - 2h[i,j]) * inv_dx2
    end
    return out
end

"""
    spatial_grid(M, N, dx, dy = dx) -> (xgrid, ygrid)

Centred coordinate vectors so that `mean(xgrid) ≈ 0`, suitable for the
wide-FOV correction `α·x·e_θ` (gauge-fixed).
"""
function spatial_grid(M::Integer, N::Integer, dx::Real, dy::Real = dx)
    xgrid = collect(((1:N) .- (N + 1) / 2) .* dx)
    ygrid = collect(((1:M) .- (M + 1) / 2) .* dy)
    return xgrid, ygrid
end
