"""
    integrate_fourier(γx, γy, κ; dx, dy = dx, ω = dx, padding = :mirror) -> δh

Step B2 of Algorithm 1: integrate the per-pixel gradient (γₓ, γᵧ) and Laplacian
(κ) into a topography correction `δh` by solving Eq. (14) of Neggers et al. 2022
in Fourier space:

    δh̃(k) = − ( i·k·γ̃(k) + ω² |k|² κ̃(k) ) / ( |k|² + ω⁴ |k|⁴ )

The DC mode is unconstrained (the gauge of `h`), so we set `δh̃(0) = 0`.

Boundary handling
- `:mirror` — reflect `γ`, `κ` to a `(2M, 2N)` grid with the symmetries that
  γₓ is odd in x / even in y, γᵧ is even in x / odd in y, and κ is even in
  both. After IFFT, return the upper-left M×N block. This kills the periodic
  wrap-around inherent in plain FFT-based integration.
- `:none` — straight FFT on the input grid (faster, but image edges leak).
"""
function integrate_fourier(γx::AbstractMatrix{T}, γy::AbstractMatrix{T},
                            κ::AbstractMatrix{T};
                            dx::Real, dy::Real = dx, ω::Real = dx,
                            padding::Symbol = :mirror) where {T}
    M, N = size(γx)
    @assert size(γy) == (M, N) && size(κ) == (M, N)
    if padding === :none
        return _integrate_fourier_core(γx, γy, κ; dx = dx, dy = dy, ω = ω)
    elseif padding === :mirror
        γx_pad = mirror_pad(γx; sym_x = :odd,  sym_y = :even)
        γy_pad = mirror_pad(γy; sym_x = :even, sym_y = :odd)
        κ_pad  = mirror_pad(κ;  sym_x = :even, sym_y = :even)
        δh_pad = _integrate_fourier_core(γx_pad, γy_pad, κ_pad; dx = dx, dy = dy, ω = ω)
        return δh_pad[1:M, 1:N]
    else
        error("padding must be :mirror or :none, got $padding")
    end
end

"""
    mirror_pad(A; sym_x, sym_y) -> A_pad

Extend `A` from `(M, N)` to `(2M, 2N)` by reflection. `sym_x` controls
the symmetry in the column direction (`:even` for f(2N+1-j) = +f(j),
`:odd` for f(2N+1-j) = -f(j)`); `sym_y` is analogous for rows.
"""
function mirror_pad(A::AbstractMatrix{T}; sym_x::Symbol, sym_y::Symbol) where {T}
    M, N = size(A)
    sx = sym_x === :even ? one(T) : sym_x === :odd ? -one(T) : error("sym_x ∈ {:even,:odd}")
    sy = sym_y === :even ? one(T) : sym_y === :odd ? -one(T) : error("sym_y ∈ {:even,:odd}")
    out = Matrix{T}(undef, 2M, 2N)
    @inbounds for j in 1:N, i in 1:M
        out[i, j]            = A[i, j]
        out[i, 2N + 1 - j]   = sx * A[i, j]
        out[2M + 1 - i, j]   = sy * A[i, j]
        out[2M + 1 - i, 2N + 1 - j] = sx * sy * A[i, j]
    end
    return out
end

function _integrate_fourier_core(γx::AbstractMatrix{T}, γy::AbstractMatrix{T},
                                  κ::AbstractMatrix{T};
                                  dx::Real, dy::Real, ω::Real) where {T}
    M, N = size(γx)
    γ̃x = fft(γx)
    γ̃y = fft(γy)
    κ̃  = fft(κ)
    kx = T.(fftfreq(N, 1 / dx) .* 2π)        # length N
    ky = T.(fftfreq(M, 1 / dy) .* 2π)        # length M
    ω² = T(ω^2); ω⁴ = T(ω^4)
    δh̃ = similar(γ̃x)
    @inbounds for j in 1:N, i in 1:M
        kxj = kx[j]; kyi = ky[i]
        k2  = kxj^2 + kyi^2
        if k2 == 0
            δh̃[i, j] = zero(eltype(δh̃))
        else
            num = im * (kxj * γ̃x[i, j] + kyi * γ̃y[i, j]) + ω² * k2 * κ̃[i, j]
            den = k2 + ω⁴ * k2^2
            δh̃[i, j] = -num / den
        end
    end
    return real.(ifft(δh̃))
end
