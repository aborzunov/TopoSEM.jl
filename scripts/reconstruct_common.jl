# Shared reconstruction pipeline used by every `scripts/reconstruct_<sample>.jl`.
#
# Each per-sample script bakes in only the fixed inputs (source TIFF directory,
# pixel size from a scale-bar measurement, calibration target, detector-angle
# convention) and delegates the actual work — load → topo_sem → detrend →
# calibrate → write CSV+TOML — to `reconstruct_and_save`.

using TopoSEM
using Statistics
using Dates

include(joinpath(@__DIR__, "csv_io.jl"))

const _PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const _DATA_DIR     = joinpath(_PROJECT_ROOT, "data")
const _DEFAULT_RECONSTRUCTION_KW = (
    # poly_order = 2 (Eq. 6 of Neggers 2022 expanded to quadratic in φ) is
    # the empirical sweet spot on our sample suite: order 1 plateaus on
    # steep-slope features (Vickers indent walls, sphere edges) leaving a
    # visibly under-fit reconstruction, while order ≥ 3 is unstable
    # without explicit ridge regularisation in step C. See README §
    # "Polynomial-order observation" for the side-by-side comparison.
    poly_order        = 2,
    use_laplacian     = false,
    use_wide_fov      = false,
    max_iter          = 60,
    tol               = 1e-4,
    weight_d          = 0.0,
    weight_σ          = 0.0,
    fourier_padding   = :mirror,
    verbose           = true,
)

"""
    reconstruct_and_save(; sample, image_files, detector_angles,
                          pixel_size_um, calibration, description = "",
                          reconstruction_kwargs = (;),
                          output_dir = data/) -> (h_μm, scale, result)

Run the full reconstruction pipeline on a 4-quadrant BSE acquisition and
persist the result to `<output_dir>/<sample>_h.csv.gz` plus
`<output_dir>/<sample>_meta.toml`.

# Arguments
- `sample` — short identifier; used as the filename prefix in `data/`.
- `image_files` — explicit list of 4 image paths in detector order
  (entry `k` corresponds to detector angle `detector_angles[k]`).
  Supports `.tif` / `.tiff` / `.jpg` / `.jpeg`.
- `detector_angles` — vector of four polar angles (rad) in load order.
- `pixel_size_um` — physical pixel size in micrometres.
- `calibration` — either `nothing` (gauge-free output) or
  `(method = :peak_to_trough | :peak_to_baseline, value = X.X)`.
- `description` — free-form text written to the meta TOML.
- `reconstruction_kwargs` — overrides over `_DEFAULT_RECONSTRUCTION_KW`,
  forwarded to `topo_sem`.
"""
function reconstruct_and_save(; sample::AbstractString,
                                image_files::AbstractVector{<:AbstractString},
                                detector_angles::AbstractVector,
                                pixel_size_um::Real,
                                calibration,
                                description::AbstractString = "",
                                reconstruction_kwargs = (;),
                                output_dir::AbstractString = _DATA_DIR)
    isdir(output_dir) || mkpath(output_dir)
    @assert length(image_files) == length(detector_angles) == 4 (
        "expected 4 image files and 4 detector angles, got " *
        "$(length(image_files)) and $(length(detector_angles))")

    println("=" ^ 70)
    println("Reconstructing sample: $sample")
    println("Files:")
    for (k, f) in enumerate(image_files)
        println("  detector $k (θ = $(round(detector_angles[k], digits=3))): $f")
    end
    println("=" ^ 70)

    images = load_bse_images(image_files)
    crop_rows = size(images[1], 1)
    println("Loaded ", length(images), " images of size ", size(images[1]),
            " (info-bar auto-cropped)")

    kw = merge(_DEFAULT_RECONSTRUCTION_KW, reconstruction_kwargs)
    result = topo_sem(images, detector_angles;
                      pixel_size = pixel_size_um,
                      ω          = pixel_size_um,
                      kw...)

    h_raw   = result.h .- (sum(result.h) / length(result.h))
    h_dt, β = detrend_planar(h_raw)

    h_out, scale = if calibration === nothing
        (h_dt, one(eltype(h_dt)))
    else
        calibrate_height(h_dt, calibration.value; method = calibration.method)
    end

    csv_path  = joinpath(output_dir, "$(sample)_h.csv.gz")
    meta_path = joinpath(output_dir, "$(sample)_meta.toml")
    write_topography_csv_gz(csv_path, h_out)

    meta = _build_meta(sample, description, collect(image_files),
                        pixel_size_um, detector_angles, crop_rows,
                        calibration, scale, β,
                        result, kw, h_out)
    write_meta_toml(meta_path, meta)

    println()
    println("Wrote topography → $(csv_path)  ($(round(filesize(csv_path)/1024, digits=1)) KB)")
    println("Wrote metadata   → $(meta_path)")
    println("h range:  [$(round(minimum(h_out), digits=4)), $(round(maximum(h_out), digits=4))]" *
            (calibration === nothing ? " (gauge-free)" : " μm"))
    println("converged = $(result.converged), iters = $(result.iterations)")

    return h_out, scale, result
end

function _build_meta(sample, description, source_files,
                     pixel_size_um, detector_angles, crop_rows,
                     calibration, scale, planar_β,
                     result, kw, h_out)
    M, N = size(h_out)
    calibration_section = if calibration === nothing
        Dict("method" => "none", "detrend_planar" => true,
             "planar_tilt_row" => planar_β[1], "planar_tilt_col" => planar_β[2])
    else
        Dict("method"          => string(calibration.method),
             "known_value_um"  => Float64(calibration.value),
             "scale_factor"    => Float64(scale),
             "detrend_planar"  => true,
             "planar_tilt_row" => planar_β[1],
             "planar_tilt_col" => planar_β[2])
    end

    final_dh = isempty(result.history) ? 0.0 : result.history[end].δh
    final_df = isempty(result.history) ? 0.0 : result.history[end].δf

    return Dict(
        "sample"        => sample,
        "description"   => description,
        "source_files"  => source_files,
        "geometry"      => Dict(
            "pixel_size_um"        => Float64(pixel_size_um),
            "detector_angles_rad"  => Float64.(collect(detector_angles)),
            "crop_rows"            => crop_rows,
        ),
        "calibration"   => calibration_section,
        "reconstruction" => Dict(
            "poly_order"      => Int(kw.poly_order),
            "use_laplacian"   => Bool(kw.use_laplacian),
            "use_wide_fov"    => Bool(kw.use_wide_fov),
            "max_iter"        => Int(kw.max_iter),
            "tol"             => Float64(kw.tol),
            "weight_sigma"    => Float64(kw.weight_σ),
            "fourier_padding" => string(kw.fourier_padding),
        ),
        "result"        => Dict(
            "size_rows"      => M,
            "size_cols"      => N,
            "h_min"          => Float64(minimum(h_out)),
            "h_max"          => Float64(maximum(h_out)),
            "h_std"          => Float64(std(vec(h_out))),
            "final_residual" => Float64(final_df),
            "final_dh_norm"  => Float64(final_dh),
            "iterations"     => result.iterations,
            "converged"      => result.converged,
        ),
        "provenance"    => Dict(
            "julia_version"   => string(VERSION),
            "toposem_version" => "0.1.0",
            "timestamp"       => string(now()),
        ),
    )
end
