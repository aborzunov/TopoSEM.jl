"""
    pixel_weights(d::DetectorModel, gx, gy; weight_d = 0.0, weight_σ = 0.5,
                  saturation_floor = 0.01) -> Matrix

Per-detector pixel weight from Eq. (15) of Neggers et al. 2022:
    w = exp(-((∂_n h − d)² + (∂_t h)²) / (2σ²))
where (∂_n, ∂_t) is the gradient rotated into the detector frame. Pixels
flagged as saturated in `d.saturated` are multiplied by `saturation_floor`
(default 1 %) so they remain in the redundant LSQ but cannot dominate it.
"""
function pixel_weights(d::DetectorModel{T}, gx::AbstractMatrix, gy::AbstractMatrix;
                       weight_d::Real = 0.0, weight_σ::Real = 0.0,
                       saturation_floor::Real = 0.01) where {T}
    M, N = size(gx)
    @assert size(gy) == (M, N)
    @assert size(d.saturated) == (M, N)
    out = Matrix{T}(undef, M, N)
    cθ = cos(d.θ); sθ = sin(d.θ)
    sf = T(saturation_floor)
    use_gaussian = weight_σ > 0
    if use_gaussian
        # Normalise slopes by the empirical scale so Eq. 15 stays meaningful
        # regardless of the absolute (gauge-free) scale of h.
        scale = sqrt(mean(gx.^2) + mean(gy.^2))
        scale = scale > eps(T) ? T(scale) : one(T)
        inv_2σ² = T(1 / (2 * (weight_σ * scale)^2))
        d_ref = T(weight_d) * scale
        @inbounds for j in 1:N, i in 1:M
            n =  cθ * gx[i,j] + sθ * gy[i,j]
            t = -sθ * gx[i,j] + cθ * gy[i,j]
            w = exp(-((n - d_ref)^2 + t^2) * inv_2σ²)
            if d.saturated[i,j]
                w *= sf
            end
            out[i,j] = w
        end
    else
        @inbounds for j in 1:N, i in 1:M
            w = one(T)
            if d.saturated[i,j]
                w *= sf
            end
            out[i,j] = w
        end
    end
    return out
end
