"""
    TopoSEMResult{T,Order}

Output of [`topo_sem`](@ref). `h` is the reconstructed topography (gauge-free
constant; subtract the mean before plotting), `detectors` holds the calibrated
response models, `residuals[k]` is the per-detector image residual at the
last iteration, `weights[k]` is the per-detector pixel weight (Eq. 15), and
`history` lists `(iter, δh, δf, ratio)` for each iteration.
"""
struct TopoSEMResult{T,Order}
    h::Matrix{T}
    detectors::Vector{DetectorModel{T,Order}}
    residuals::Vector{Matrix{T}}
    weights::Vector{Matrix{T}}
    history::Vector{NamedTuple{(:iter, :δh, :δf, :ratio), NTuple{4,T}}}
    converged::Bool
    iterations::Int
end

"""
    topo_sem(images, angles; kwargs...) -> TopoSEMResult

Reconstruct surface topography `h(x, y)` from `N ≥ 3` co-registered BSE
detector images (Algorithm 1 of Neggers et al. 2022).

# Arguments
- `images::Vector{<:AbstractMatrix}` — `N` images of the same size, in
  arbitrary intensity units. The conventional 0..255 (8-bit) raw range is
  fine; pass `bits = 8` so saturation is detected at 0 and 255.
- `angles::AbstractVector` — `N` polar angles `θ⁽ⁱ⁾` of the detector centres
  of mass in the image frame (radians). For a 4-quadrant detector the typical
  convention is `(0, π/2, π, 3π/2)`.

# Keyword arguments
- `pixel_size = 1.0` — physical pixel size (any consistent unit).
- `ω = pixel_size` — gradient/curvature crossover wavelength in Eq. 14.
- `poly_order = 1` — total degree of the bivariate polynomial `φ` (≥ 1).
- `use_laplacian = false` — include the `p₃·∇²h` term.
- `use_wide_fov = false` — include the wide-FOV correction `c(x, h)` (α, β).
- `max_iter = 30`, `tol = 1e-4` — outer-loop budget and relative `‖δh‖/‖h‖`
  convergence threshold.
- `weight_d = 0.0`, `weight_σ = 0.5` — Eq. 15 reference normal slope and
  confidence width. Defaults give a near-uniform weight on flat regions.
- `bits = 8` — bit depth of the input images, used for saturation detection.
- `verbose = true` — print per-iteration log to stdout.
- `fourier_padding = :mirror` — `:mirror` or `:none` (see [`integrate_fourier`](@ref)).
"""
function topo_sem(images::Vector{<:AbstractMatrix}, angles::AbstractVector;
                  pixel_size::Real = 1.0,
                  ω::Real = pixel_size,
                  poly_order::Int = 1,
                  use_laplacian::Bool = false,
                  use_wide_fov::Bool = false,
                  max_iter::Int = 30,
                  tol::Real = 1e-4,
                  weight_d::Real = 0.0,
                  weight_σ::Real = 0.0,
                  bits::Int = 8,
                  verbose::Bool = true,
                  fourier_padding::Symbol = :mirror)

    @assert length(images) == length(angles) "`images` and `angles` must match in length"
    @assert length(images) ≥ 3 "Algorithm requires N ≥ 3 detectors"
    T = Float64
    M, N = size(images[1])
    for k in eachindex(images)
        @assert size(images[k]) == (M, N) "image $k has size $(size(images[k])); expected $((M, N))"
    end

    Order = poly_order
    max_val = T(2^bits - 1)

    img_T = [T.(im) for im in images]
    detectors = [DetectorModel(T(angles[k]), img_T[k]; order = Order,
                               max_val = max_val, eltype = T)
                 for k in eachindex(img_T)]

    h    = zeros(T, M, N)
    lap  = zeros(T, M, N)
    γx   = zeros(T, M, N)
    γy   = zeros(T, M, N)
    κ    = zeros(T, M, N)
    xgrid, ygrid = spatial_grid(M, N, T(pixel_size))

    residuals = [Matrix{T}(undef, M, N) for _ in 1:length(detectors)]
    weights   = [Matrix{T}(undef, M, N) for _ in 1:length(detectors)]
    pred      = Matrix{T}(undef, M, N)

    history = NamedTuple{(:iter, :δh, :δf, :ratio), NTuple{4,T}}[]
    converged = false
    final_iter = 0

    dx  = T(pixel_size)
    ω_T = T(ω)

    for iter in 1:max_iter
        final_iter = iter

        gx, gy = gradient_centered(h, dx)
        if use_laplacian
            lap = laplacian_5pt(h, dx)
        end

        for k in eachindex(detectors)
            weights[k] .= pixel_weights(detectors[k], gx, gy;
                                        weight_d = weight_d, weight_σ = weight_σ)
        end

        if iter > 1
            recalibrate!(detectors, img_T, h, gx, gy, lap, xgrid, ygrid, weights;
                         use_laplacian = use_laplacian,
                         use_wide_fov = use_wide_fov)
            fix_gauge!(detectors, h)
            # Re-evaluate derivatives after h rescaling
            gx, gy = gradient_centered(h, dx)
            if use_laplacian
                lap = laplacian_5pt(h, dx)
            end
        end

        δf_norm² = zero(T)
        for k in eachindex(detectors)
            forward_image!(pred, detectors[k], h, gx, gy, lap, xgrid, ygrid;
                            use_laplacian = use_laplacian,
                            use_wide_fov = use_wide_fov)
            residuals[k] .= img_T[k] .- pred
            δf_norm² += sum(abs2, residuals[k])
        end
        δf_norm = sqrt(δf_norm²)

        solve_pixelwise!(γx, γy, κ, residuals, detectors, h, gx, gy, weights;
                         use_laplacian = use_laplacian,
                         use_wide_fov = use_wide_fov,
                         xgrid = xgrid, ygrid = ygrid)

        δh = integrate_fourier(γx, γy, κ; dx = dx, dy = dx, ω = ω_T,
                                padding = fourier_padding)

        h .+= δh

        δh_norm = sqrt(sum(abs2, δh))
        h_norm  = sqrt(sum(abs2, h)) + eps(T)
        ratio   = δh_norm / h_norm
        push!(history, (iter = T(iter), δh = δh_norm, δf = δf_norm, ratio = ratio))

        if verbose
            @printf("[TopoSEM] iter %2d: ‖δh‖ = %.4e   ‖δf‖ = %.4e   ratio = %.3e\n",
                    iter, δh_norm, δf_norm, ratio)
        end

        if iter ≥ 2 && ratio < T(tol)
            converged = true
            break
        end
    end

    return TopoSEMResult{T,Order}(h, detectors, residuals, weights,
                                   history, converged, final_iter)
end

"""
    calibrate_height(h, known_height_μm; method = :peak_to_baseline) -> (h_μm, scale)

Convert a gauge-free reconstructed topography `h` to physical micrometres by
matching one feature of known dimension. Without this step the reconstructed
heights are dimensionless (only relative shapes are meaningful), since the
forward model is invariant under `(p₂ → p₂/c, h → h·c)`.

`method`:
- `:peak_to_baseline` — assume `known_height_μm` is the difference between
  the highest reconstructed point (the apex) and the median (the surrounding
  reference plane). Robust to small saturation / boundary artefacts.
- `:peak_to_trough`   — assume `known_height_μm = max(h) − min(h)`. Use
  this if both extremes are real features.

Returns `(h_μm, scale)` where `h_μm = h .* scale` and `scale` is the
multiplicative factor applied. For the Pyramid sample in Neggers et al. 2022,
`known_height_μm ≈ 2.0` (three steps totalling 600 + 600 + 800 nm). For the
Dome / sphere sample, `known_height_μm ≈ 1.5`.
"""
function calibrate_height(h::AbstractMatrix, known_height_μm::Real;
                          method::Symbol = :peak_to_baseline)
    if method === :none
        return h, one(eltype(h))
    elseif method === :peak_to_baseline
        baseline = median(h)
        observed = maximum(h) - baseline
    elseif method === :peak_to_trough
        observed = maximum(h) - minimum(h)
    else
        error("method must be :none, :peak_to_baseline or :peak_to_trough; got $method")
    end
    @assert observed > 0 "reconstructed h is degenerate (max = $(maximum(h)), baseline = $(method === :peak_to_baseline ? median(h) : minimum(h)))"
    scale = known_height_μm / observed
    return h .* scale, scale
end

"""
    detrend_planar(h) -> (h_detrend, β)

Subtract a least-squares planar tilt `h ≈ a·i + b·j + c` (with `i`, `j` row
and column indices) from `h`. Returns the detrended array and the fit
coefficients `β = [a, b, c]`. Useful before visualisation when stage tilt
or detector-bias asymmetry dominates over the topographic signal of
interest (small features on top of an essentially flat surface).
"""
function detrend_planar(h::AbstractMatrix{T}) where {T}
    M, N = size(h)
    ii = repeat(collect(T(1):T(M)), 1, N)
    jj = repeat(collect(T(1):T(N))', M, 1)
    A  = hcat(vec(ii), vec(jj), ones(T, M*N))
    β  = A \ vec(h)
    h_dt = h .- reshape(A * β, M, N)
    return h_dt, β
end
