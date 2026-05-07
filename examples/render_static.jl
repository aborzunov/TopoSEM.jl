# Render a static PNG figure (heatmap + 3D surface) from a pre-computed
# reconstruction. Used by `make figures`.
#
# Usage:
#     julia --project=. examples/render_static.jl <sample>
#
# Output: `data/<sample>_topography.png`.

using CairoMakie

include(joinpath(@__DIR__, "..", "scripts", "csv_io.jl"))

const _PROJECT_ROOT = normpath(joinpath(@__DIR__, ".."))
const _DATA_DIR     = joinpath(_PROJECT_ROOT, "data")
const _SAMPLES = ("jun_vickers", "jun_vickers_exp1",
                   "jun_sphere", "jun_sphere_exp2", "jun_sphere_exp3",
                   "feb_P", "feb_PR", "feb_S", "feb_V",
                   "may_sphere", "may_piezoceramic",
                   "may_crooked_60", "may_crooked_30")

function _parse_sample(args)
    if isempty(args)
        println(stderr, "usage: julia examples/render_static.jl <sample>")
        exit(1)
    end
    sample = args[1]
    sample ∈ _SAMPLES || error("unknown sample '$sample'; expected one of $(_SAMPLES)")
    return sample
end

function render(sample::AbstractString)
    csv_path  = joinpath(_DATA_DIR, "$(sample)_h.csv.gz")
    meta_path = joinpath(_DATA_DIR, "$(sample)_meta.toml")
    isfile(csv_path) || error("missing $(csv_path) — run `make $(sample)` first")

    h = read_topography_csv_gz(csv_path)
    meta = read_meta_toml(meta_path)
    M, N = size(h)

    pixel_size = meta["geometry"]["pixel_size_um"]
    xs = ((1:N) .- (N + 1) / 2) .* pixel_size
    ys = ((1:M) .- (M + 1) / 2) .* pixel_size
    Lx = xs[end] - xs[1]; Ly = ys[end] - ys[1]

    cal = meta["calibration"]
    z_label = cal["method"] == "none" ? "h (gauge-free)" : "h (μm)"
    title_calib = cal["method"] == "none" ? "uncalibrated" :
                  "$(cal["method"]) = $(cal["known_value_um"]) μm"

    fig = Figure(; size = (1500, 600))
    ax1 = Axis(fig[1, 1]; aspect = DataAspect(),
               xlabel = "x (μm)", ylabel = "y (μm)",
               title  = "$(sample) — top-down ($(title_calib))")
    hm = heatmap!(ax1, xs, ys, h'; colormap = :viridis)
    Colorbar(fig[1, 2], hm, label = z_label)

    ax2 = Axis3(fig[1, 3]; aspect = (1, Ly / Lx, 0.45),
                xlabel = "x (μm)", ylabel = "y (μm)", zlabel = z_label,
                title  = "$(sample) — 3D surface")
    surface!(ax2, xs, ys, h'; colormap = :viridis)

    out_path = joinpath(_DATA_DIR, "$(sample)_topography.png")
    save(out_path, fig)
    println("Saved → $(out_path)")
end

if abspath(PROGRAM_FILE) == @__FILE__
    sample = _parse_sample(ARGS)
    render(sample)
end
