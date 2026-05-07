"""
    load_bse_image(path; bits = 8, crop = :auto) -> Matrix{Float64}

Load a single grayscale TIFF and return it as a `Matrix{Float64}` rescaled
to the raw 0..(2^bits − 1) range. This matches the convention used by the
saturation detection in [`DetectorModel`](@ref) and by the default `bits = 8`
keyword on [`topo_sem`](@ref).

`crop` controls the removal of the SEM info bar (the dark plate at the
bottom of the image with magnification, dwell time, etc. — non-physical
intensity that corrupts the LSQ residual):
- `:auto` (default) — detect the plate from the sharp drop in row-wise
  intensity standard deviation and crop it.
- `:none` or `nothing` — return the full image as stored in the TIFF.
- `UnitRange{Int}` — explicit vertical row range to keep, e.g. `1:718`.
"""
function load_bse_image(path::AbstractString; bits::Int = 8,
                        crop::Union{Symbol,Nothing,UnitRange{<:Integer}} = :auto)
    img = _load_image_native(path)
    scale = float(2^bits - 1)
    out = Matrix{Float64}(undef, size(img)...)
    @inbounds for i in eachindex(img)
        out[i] = Float64(img[i]) * scale
    end
    rows = _resolve_crop(out, crop)
    return out[rows, :]
end

"""
    _load_image_native(path) -> Matrix{<:Gray}

Dispatch on file extension: TIFF via `TiffImages.jl`, JPEG via `FileIO + JpegTurbo`.
RGB JPEGs (some SEMs save colour artifacts of the BSE quadrant overlay) are
converted to grayscale with a luminance-weighted mean.
"""
function _load_image_native(path::AbstractString)
    ext = lowercase(splitext(path)[2])
    if ext in (".tif", ".tiff")
        return TiffImages.load(path)
    elseif ext in (".jpg", ".jpeg")
        img = FileIO.load(path)
        return _to_gray(img)
    else
        error("unsupported image extension '$(ext)' for '$(path)'; expected .tif/.tiff/.jpg/.jpeg")
    end
end

# Already grayscale — pass through.
_to_gray(img::AbstractMatrix{<:Real}) = img
_to_gray(img::AbstractMatrix{<:Gray}) = img
# RGB → grayscale with Rec. 709 luminance weights.
function _to_gray(img::AbstractMatrix)
    eltype_str = string(eltype(img))
    occursin("RGB", eltype_str) || error(
        "unexpected pixel type $(eltype(img)) — expected Gray or RGB")
    return Gray.(img)
end

_resolve_crop(::AbstractMatrix, ::Nothing)    = (:)
_resolve_crop(::AbstractMatrix, c::UnitRange) = c
function _resolve_crop(img::AbstractMatrix, sym::Symbol)
    if sym === :none
        return 1:size(img, 1)
    elseif sym === :auto
        return detect_data_rows(img)
    else
        error("crop must be :auto, :none, nothing, or UnitRange; got $sym")
    end
end

"""
    detect_data_rows(img; mean_dev_threshold = 50.0, std_threshold_factor = 0.5,
                     search_fraction = 0.4, min_data_rows = 200,
                     safety_margin = 5)

Heuristically find pure-data rows at the top of an SEM acquisition by
detecting any row in the bottom `search_fraction` whose statistics deviate
from the data baseline. SEM overlays the bottom of an image with a
heterogeneous mix of artefacts (any of which can sit *inside* what would
otherwise be pure data):
  - dark divider lines (mean ≈ 0, very low std)
  - bright scale bars (mean ≈ saturation, low std)
  - text plates (mean = arbitrary, low std)
  - mixed bands (e.g. text plus bar)

A row is flagged as non-data if **either**:
  - `|mean(row) − ref_mean|` exceeds `mean_dev_threshold` in raw 0..255
    grayscale units (catches saturated / blacked-out overlays);
  - `std(row) < std_threshold_factor × ref_std` (catches uniform plates).

Reference (`ref_mean`, `ref_std`) is taken from the top half of the image,
which is assumed to be pure data. The kept range is
`1:(top_non_data − 1 − safety_margin)` so even the highest detected
anomaly in the bottom `search_fraction` is excluded.
"""
function detect_data_rows(img::AbstractMatrix; mean_dev_threshold::Real = 50.0,
                           std_threshold_factor::Real = 0.5,
                           search_fraction::Real = 0.15,
                           min_data_rows::Int = 200,
                           safety_margin::Int = 5)
    M = size(img, 1)
    row_mean = [mean(view(img, i, :)) for i in 1:M]
    row_std  = [std(view(img, i, :))  for i in 1:M]
    top_n = max(1, M ÷ 2)
    ref_mean = median(view(row_mean, 1:top_n))
    ref_std  = median(view(row_std,  1:top_n))
    std_threshold = std_threshold_factor * max(ref_std, eps(eltype(row_std)))
    search_start = max(min_data_rows + 1, Int(round((1 - search_fraction) * M)))
    top_non_data = M + 1
    for i in M:-1:search_start
        is_anomalous = abs(row_mean[i] - ref_mean) > mean_dev_threshold ||
                       row_std[i] < std_threshold
        if is_anomalous
            top_non_data = i
        end
    end
    last_data = top_non_data == M + 1 ? M :
                max(min_data_rows, top_non_data - 1 - safety_margin)
    return 1:last_data
end

"""
    load_bse_quadrants(dir; pattern = r"^[QS]([1-4])\\.tif\$",
                       angles = (0.0, π/2, π, 3π/2),
                       quadrant_permutation = (1, 2, 3, 4),
                       flip_y = false,
                       bits = 8) -> (Vector{Matrix{Float64}}, Vector{Float64})

Load four BSE detector images named `Q1..Q4` (or `S1..S4` for the sphere
sample) from `dir` and return them along with their assumed polar angles.

The default `angles` correspond to the convention in Eq. 10 of Neggers et al.
2022: `f⁽¹⁾` is sensitive to `+∂ₓh`, `f⁽²⁾` to `+∂ᵧh`, etc. If the actual SEM
quadrant numbering differs, override `quadrant_permutation` (a permutation of
`(1,2,3,4)` saying which input file should be treated as Q1, Q2, Q3, Q4) or
`flip_y = true` (mirror the image rows when the SEM y-axis is inverted).
"""
function load_bse_quadrants(dir::AbstractString;
                            pattern::Regex = r"^[QS]([1-4])\.tif$",
                            angles = (0.0, π/2, π, 3π/2),
                            quadrant_permutation::NTuple{4,Int} = (1, 2, 3, 4),
                            flip_y::Bool = false,
                            bits::Int = 8,
                            crop::Union{Symbol,Nothing,UnitRange{<:Integer}} = :auto)
    @assert isdir(dir) "directory not found: $dir"
    @assert sort(collect(quadrant_permutation)) == [1, 2, 3, 4] "permutation must be of (1,2,3,4)"

    files = sort(readdir(dir))
    # Sort by filename — works for both numeric (Q1..Q4) and alphabetic
    # (PA..PD) quadrant suffixes, since both lexicographically order the way
    # we want. For arbitrary numeric ordering (e.g. May series with files
    # 13/16/1/4), use [`load_bse_images`](@ref) with explicit paths.
    matched_files = sort([f for f in files if occursin(pattern, f)])
    @assert length(matched_files) == 4 "expected 4 files matching $pattern in $dir, got $(length(matched_files))"
    matched_paths = [joinpath(dir, f) for f in matched_files]
    raw_full = [load_bse_image(p; bits = bits, crop = :none) for p in matched_paths]
    # For :auto, take the MOST CONSERVATIVE (smallest last_data) crop across
    # all four quadrants — overlay artefacts may extend higher in some
    # detectors than others, and we want all four images to share the same
    # row range so the LSQ inputs are co-registered.
    if crop === :auto
        last_data = minimum(last(detect_data_rows(im)) for im in raw_full)
        rows = 1:last_data
    else
        rows = _resolve_crop(raw_full[1], crop)
    end
    raw = [im[rows, :] for im in raw_full]
    if flip_y
        raw = [reverse(im, dims = 1) for im in raw]
    end
    images = [raw[quadrant_permutation[k]] for k in 1:4]
    return images, collect(Float64, angles)
end

"""
    load_bse_images(paths; bits = 8, crop = :auto) -> Vector{Matrix{Float64}}

Load N BSE images from an **explicit list of paths**, in the order given,
and apply a co-registered crop (most-conservative `last_data` across all
channels for `:auto`). Use this when filenames don't sort into the desired
detector order — e.g., the May series where the four quadrants of one
acquisition are arbitrary indices like `13.jpg, 16.jpg, 1.jpg, 4.jpg`.

For lexicographically-ordered groups (Q1..Q4, PA..PD) prefer
[`load_bse_quadrants`](@ref) which reads from a directory + regex.
"""
function load_bse_images(paths::AbstractVector{<:AbstractString};
                          bits::Int = 8,
                          crop::Union{Symbol,Nothing,UnitRange{<:Integer}} = :auto)
    @assert !isempty(paths) "at least one image path required"
    raw_full = [load_bse_image(p; bits = bits, crop = :none) for p in paths]
    if crop === :auto
        last_data = minimum(last(detect_data_rows(im)) for im in raw_full)
        rows = 1:last_data
    else
        rows = _resolve_crop(raw_full[1], crop)
    end
    return [im[rows, :] for im in raw_full]
end
