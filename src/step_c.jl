"""
    recalibrate!(detectors, images, h, gx, gy, lap, xgrid, ygrid, weights;
                 use_laplacian = false, use_wide_fov = false,
                 inner_iter = 3) -> detectors

Step C of Algorithm 1: with `h` (and its derivatives) held fixed, refit each
detector's response model by weighted least squares. The free parameters per
detector are

    p₀, p₂ₙ, p₂ₜ, [p₃,]  [A, B,]  q_ij…

with `A = p₂ₙ·α`, `B = p₂ₙ·β` (re-parametrisation that keeps the LSQ linear
even with the wide-FOV coupling — `α, β` are recovered as `A/p₂ₙ`, `B/p₂ₙ`).

The polynomial monomials use `arg_n_prev` (i.e., the current `(α, β)` when
`use_wide_fov`); for `Order ≥ 2` we therefore loop `inner_iter` times so
`(α, β)` and the higher-order coefficients converge together.
"""
function recalibrate!(detectors::Vector{<:DetectorModel{T}},
                      images::Vector{<:AbstractMatrix},
                      h::AbstractMatrix, gx::AbstractMatrix, gy::AbstractMatrix,
                      lap::AbstractMatrix,
                      xgrid::AbstractVector, ygrid::AbstractVector,
                      weights::Vector{<:AbstractMatrix};
                      use_laplacian::Bool = false,
                      use_wide_fov::Bool = false,
                      inner_iter::Int = 3) where {T}
    M, N = size(h)
    K_lin = 3 + (use_laplacian ? 1 : 0) + (use_wide_fov ? 2 : 0)

    for k in eachindex(detectors)
        d = detectors[k]
        Order = detector_order(d)
        nh = n_higher_coeffs(Order)
        K = K_lin + nh

        AtWA = zeros(T, K, K)
        AtWr = zeros(T, K)
        row  = zeros(T, K)
        mono = nh > 0 ? zeros(T, nh) : T[]

        # Sub-iterations matter only when the design matrix entries depend on
        # the parameters being fit, which is the case for Order ≥ 2 with
        # wide-FOV (the monomials of `arg_n` then drift with α, β).
        n_inner = (Order ≥ 2 && use_wide_fov) ? inner_iter : 1

        for _ in 1:n_inner
            fill!(AtWA, zero(T)); fill!(AtWr, zero(T))
            cθ = cos(d.θ); sθ = sin(d.θ)
            @inbounds for j in 1:N, i in 1:M
                G   =  cθ * gx[i,j] + sθ * gy[i,j]
                G_t = -sθ * gx[i,j] + cθ * gy[i,j]
                X_e = xgrid[j] * cθ + ygrid[i] * sθ
                # Argument used to evaluate higher-order monomials (stale α, β).
                n_arg = G
                if use_wide_fov
                    n_arg += d.α * X_e + d.β * h[i,j]
                end
                t_arg = G_t

                # Design row
                row[1] = one(T)            # p₀
                row[2] = G                  # p₂ₙ
                row[3] = G_t                # p₂ₜ
                idx = 4
                if use_laplacian
                    row[idx] = lap[i,j]; idx += 1   # p₃
                end
                if use_wide_fov
                    row[idx] = X_e;       idx += 1   # A = p₂ₙ·α
                    row[idx] = h[i,j];    idx += 1   # B = p₂ₙ·β
                end
                if nh > 0
                    monomial_basis_higher!(mono, d.poly, n_arg, t_arg)
                    for kk in 1:nh
                        row[idx] = mono[kk]; idx += 1
                    end
                end

                w = weights[k][i,j]
                r = images[k][i,j]
                for ii in 1:K
                    AtWr[ii] += w * row[ii] * r
                    @simd for jj in 1:K
                        AtWA[ii,jj] += w * row[ii] * row[jj]
                    end
                end
            end

            # Mild Tikhonov regularisation against rank deficiency on flat regions.
            ε = T(1e-10) * (sum(AtWA[ii,ii] for ii in 1:K) / K + eps(T))
            for ii in 1:K
                AtWA[ii,ii] += ε
            end

            sol = AtWA \ AtWr

            d.p₀ = sol[1]
            d.p₂ = SVector(sol[2], sol[3])
            idx = 4
            if use_laplacian
                d.p₃ = sol[idx]; idx += 1
            end
            if use_wide_fov
                A_par = sol[idx]; idx += 1
                B_par = sol[idx]; idx += 1
                p2n = d.p₂[1]
                if abs(p2n) > T(1e-12)
                    d.α = A_par / p2n
                    d.β = B_par / p2n
                end
            end
            for kh in 1:nh
                d.poly.coeffs[kh] = sol[idx]; idx += 1
            end
        end
    end
    return detectors
end

"""
    fix_gauge!(detectors, h) -> s

Resolve the multiplicative gauge of the response model: rescale all detectors
so that `mean(norm(p₂))_i = 1` and rescale `h` (in place) by the same factor
so that the forward model `p₂·∇h` is invariant. The other parameters scale
in the dimensionally-consistent way:
    p₂ ← p₂/s,   p₃ ← p₃/s,   α ← α·s,   β ← β,
    higher-order polynomial: `poly.scale ← poly.scale·s`,
    h  ← h·s.
Without this step the iteration has a free dimension along which it can
drift indefinitely (manifests as oscillating ‖δh‖).
"""
function fix_gauge!(detectors::Vector{<:DetectorModel{T,Order}},
                    h::AbstractMatrix) where {T,Order}
    norms = T[norm(d.p₂) for d in detectors]
    s = mean(norms)
    if !(s > T(1e-12)) || !isfinite(s)
        return one(T)
    end
    for d in detectors
        d.p₂ = d.p₂ / s
        d.p₃ = d.p₃ / s
        d.α  = d.α  * s
        if Order ≥ 2
            d.poly = BivariatePoly{Order,T}(copy(d.poly.coeffs), d.poly.scale * s)
        end
    end
    h .*= s
    return s
end
