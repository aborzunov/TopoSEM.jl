"""
    DetectorModel{T,Order}

Per-detector parameters of the BSE response model from Neggers et al. 2022.

Conventions
- `θ` is the initial polar angle of the detector centre of mass in the world
  frame (image axes). It defines a rotation `R(θ)` that maps world gradients
  `(∂ₓh, ∂ᵧh)` into the detector-aligned frame `(n, t)` where `n = e_θ`.
- `p₂ = SVector(p₂_n, p₂_t)` lives in the detector frame. Initialisation is
  `(1, 0)` (Algorithm 1 sets `p₂ = e_θ` in the world frame, which becomes
  `(1, 0)` after rotation by `R(θ)`).
- `p₃` multiplies the Laplacian `∇²h`.
- `α`, `β` are the wide-field-of-view corrections from Eq. 4–5; the argument
  of `φ` becomes `arg_n = (R(θ)·∇h)_n + α·(x·e_θ) + β·h`, `arg_t = (R(θ)·∇h)_t`.
- `poly` holds the higher-order (degree ≥ 2) part of φ.
- `saturated` masks pixels at the extremes of the dynamic range (low weight).
"""
mutable struct DetectorModel{T,Order}
    θ::T
    p₀::T
    p₂::SVector{2,T}
    p₃::T
    α::T
    β::T
    poly::BivariatePoly{Order,T}
    saturated::BitMatrix
end

function DetectorModel(θ::Real, image::AbstractMatrix; order::Int = 1,
                        max_val::Real = 1.0, eltype::Type{T} = Float64) where {T}
    img = T.(image)
    p₀ = T(mean(img))
    p₂ = SVector(one(T), zero(T))            # detector-frame init
    p₃ = zero(T)
    α  = zero(T); β = zero(T)
    poly = BivariatePoly(T, order)
    sat_lo = img .≤ zero(T)
    sat_hi = img .≥ T(max_val)
    saturated = BitMatrix(sat_lo .| sat_hi)
    return DetectorModel{T,order}(T(θ), p₀, p₂, p₃, α, β, poly, saturated)
end

@inline detector_order(::DetectorModel{T,Order}) where {T,Order} = Order
@inline detector_eltype(::DetectorModel{T,Order}) where {T,Order} = T

"""
    rotate_to_detector(d, gx, gy) -> (arg_n, arg_t)

Rotate a world-frame gradient `(gx, gy)` into the detector frame using `θ`.
"""
@inline function rotate_to_detector(d::DetectorModel, gx::Real, gy::Real)
    c = cos(d.θ); s = sin(d.θ)
    return (c * gx + s * gy, -s * gx + c * gy)
end

"""
    detector_argument(d, gx, gy, h, x, y; use_wide_fov)

Return `(arg_n, arg_t)` — the argument of φ at one pixel.
"""
@inline function detector_argument(d::DetectorModel, gx::Real, gy::Real,
                                   h::Real, x::Real, y::Real;
                                   use_wide_fov::Bool = false)
    n, t = rotate_to_detector(d, gx, gy)
    if use_wide_fov
        c = cos(d.θ); s = sin(d.θ)
        n += d.α * (x * c + y * s) + d.β * h
    end
    return (n, t)
end

"""
    forward_image!(out, d, h, gx, gy, lap, xgrid, ygrid;
                   use_laplacian = false, use_wide_fov = false)

Evaluate the forward model `f̂^(i)(h)` (Eq. 8) over the whole grid, writing
into `out`.
"""
function forward_image!(out::AbstractMatrix, d::DetectorModel{T,Order},
                        h::AbstractMatrix, gx::AbstractMatrix, gy::AbstractMatrix,
                        lap::AbstractMatrix,
                        xgrid::AbstractVector, ygrid::AbstractVector;
                        use_laplacian::Bool = false,
                        use_wide_fov::Bool = false) where {T,Order}
    M, N = size(h)
    @assert size(out) == (M, N)
    p₂n = d.p₂[1]; p₂t = d.p₂[2]
    cθ = cos(d.θ); sθ = sin(d.θ)
    @inbounds for j in 1:N
        x = xgrid[j]
        for i in 1:M
            y = ygrid[i]
            n_arg = cθ * gx[i,j] + sθ * gy[i,j]
            t_arg = -sθ * gx[i,j] + cθ * gy[i,j]
            if use_wide_fov
                n_arg += d.α * (x * cθ + y * sθ) + d.β * h[i,j]
            end
            v = d.p₀ + p₂n * n_arg + p₂t * t_arg + evaluate_higher(d.poly, n_arg, t_arg)
            if use_laplacian
                v += d.p₃ * lap[i,j]
            end
            out[i,j] = v
        end
    end
    return out
end

function forward_image(d::DetectorModel, h, gx, gy, lap, xgrid, ygrid; kwargs...)
    out = similar(h)
    forward_image!(out, d, h, gx, gy, lap, xgrid, ygrid; kwargs...)
    return out
end

# Convenience alias matching the exported name in the plan
const forward_model = forward_image
