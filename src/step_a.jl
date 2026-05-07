"""
    solve_pixelwise!(γx, γy, κ, residuals, detectors, h, gx, gy, weights;
                     use_laplacian = false, use_wide_fov = false,
                     xgrid = nothing, ygrid = nothing) -> (γx, γy, κ)

Step A2 of Algorithm 1 (Neggers et al. 2022). At every pixel build the linear
system that comes from the tangent model (Eq. 9):

    δf̂^(i) ≈ a^(i)·γₓ + b^(i)·γᵧ + c^(i)·κ

where the coefficients depend on the local detector Jacobian
`(Jₙ, Jₜ) = ∇φ^(i)(R(θ)·∇h + c)` and the rotation angle θ:

    a^(i) =  Jₙ·cos θ − Jₜ·sin θ
    b^(i) =  Jₙ·sin θ + Jₜ·cos θ
    c^(i) =  p₃^(i)                      (only when `use_laplacian`)

The per-pixel solve is a weighted least squares of size `N × M` with `N`
the number of detectors and `M ∈ {2, 3}` unknowns (2 if `use_laplacian` is
disabled). `weights[k][i,j]` is the weight of pixel (i,j) for detector k.
"""
function solve_pixelwise!(γx::AbstractMatrix{T}, γy::AbstractMatrix{T},
                          κ::AbstractMatrix{T},
                          residuals::Vector{<:AbstractMatrix},
                          detectors::Vector{<:DetectorModel},
                          h::AbstractMatrix, gx::AbstractMatrix, gy::AbstractMatrix,
                          weights::Vector{<:AbstractMatrix};
                          use_laplacian::Bool = false,
                          use_wide_fov::Bool = false,
                          xgrid::Union{Nothing,AbstractVector} = nothing,
                          ygrid::Union{Nothing,AbstractVector} = nothing) where {T}
    M, N = size(h)
    K = length(detectors)
    @assert length(residuals) == K
    @assert length(weights) == K
    nu = use_laplacian ? 3 : 2
    if use_wide_fov
        @assert xgrid !== nothing && ygrid !== nothing
    end

    # Pre-compute trigonometric coefficients per detector
    cθs = [cos(d.θ) for d in detectors]
    sθs = [sin(d.θ) for d in detectors]

    Threads.@threads for i in 1:M
        AtWA = zeros(T, nu, nu)
        AtWr = zeros(T, nu)
        row  = zeros(T, nu)
        @inbounds for j in 1:N
            fill!(AtWA, zero(T)); fill!(AtWr, zero(T))
            for k in 1:K
                d = detectors[k]
                cθ = cθs[k]; sθ = sθs[k]
                # Argument of φ at this pixel for detector k
                n_arg =  cθ * gx[i,j] + sθ * gy[i,j]
                t_arg = -sθ * gx[i,j] + cθ * gy[i,j]
                if use_wide_fov
                    n_arg += d.α * (xgrid[j] * cθ + ygrid[i] * sθ) + d.β * h[i,j]
                end
                # Jacobian (∂φ/∂u, ∂φ/∂v)
                Jh_u, Jh_v = gradient_higher(d.poly, n_arg, t_arg)
                Jn = d.p₂[1] + Jh_u
                Jt = d.p₂[2] + Jh_v
                a = Jn * cθ - Jt * sθ
                b = Jn * sθ + Jt * cθ
                row[1] = a; row[2] = b
                if use_laplacian
                    row[3] = d.p₃
                end
                w = weights[k][i,j]
                r = residuals[k][i,j]
                # Accumulate AᵀWA and AᵀWr
                @simd for ii in 1:nu
                    AtWr[ii] += w * row[ii] * r
                    for jj in 1:nu
                        AtWA[ii,jj] += w * row[ii] * row[jj]
                    end
                end
            end
            # Per-pixel solve. Tikhonov shim avoids singularity on flat regions.
            ε = T(1e-10) * (AtWA[1,1] + AtWA[2,2] + (use_laplacian ? AtWA[3,3] : zero(T)))
            for ii in 1:nu
                AtWA[ii,ii] += ε
            end
            sol = AtWA \ AtWr
            γx[i,j] = sol[1]
            γy[i,j] = sol[2]
            κ[i,j]  = use_laplacian ? sol[3] : zero(T)
        end
    end
    return γx, γy, κ
end
