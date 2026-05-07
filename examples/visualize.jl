# Interactive 3D viewer for a pre-computed reconstruction.
#
# Usage:
#     julia --project=. examples/visualize.jl <sample>
#
# where <sample> is one of:
#     pyramid sphere feb_P feb_PR feb_S feb_V
#
# Loads `data/<sample>_h.csv.gz` and `data/<sample>_meta.toml`, then opens a
# native OpenGL window (GLMakie). Drag to rotate, scroll to zoom, right-drag
# to pan. Close the window to exit.
#
# The reconstruction is NOT re-run — only loaded. To regenerate the CSV from
# raw TIFFs, use `make <sample>` (or `make all`).

using GLMakie

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
        println(stderr, "usage: julia examples/visualize.jl <sample>")
        println(stderr, "       <sample> ∈ {", join(_SAMPLES, ", "), "}")
        exit(1)
    end
    sample = args[1]
    sample ∈ _SAMPLES || error("unknown sample '$sample'; expected one of $(_SAMPLES)")
    return sample
end

function visualize(sample::AbstractString)
    csv_path  = joinpath(_DATA_DIR, "$(sample)_h.csv.gz")
    meta_path = joinpath(_DATA_DIR, "$(sample)_meta.toml")
    isfile(csv_path)  || error("missing $(csv_path) — run `make $(sample)` first")
    isfile(meta_path) || error("missing $(meta_path) — run `make $(sample)` first")

    h = read_topography_csv_gz(csv_path)
    meta = read_meta_toml(meta_path)
    M, N = size(h)

    pixel_size = meta["geometry"]["pixel_size_um"]
    xs = ((1:N) .- (N + 1) / 2) .* pixel_size
    ys = ((1:M) .- (M + 1) / 2) .* pixel_size
    Lx = xs[end] - xs[1]; Ly = ys[end] - ys[1]

    cal = meta["calibration"]
    z_label, sub_title = if cal["method"] == "none"
        ("h (gauge-free units)",
         "uncalibrated — vertical scale arbitrary")
    else
        ("h (μm)",
         "calibrated: $(cal["method"]) = $(cal["known_value_um"]) μm")
    end

    GLMakie.activate!(; title = "TopoSEM — $(sample) (drag to rotate)")
    fig = Figure(; size = (1300, 850))
    ax = Axis3(fig[1, 1];
               aspect = (1, Ly / Lx, 0.45),
               xlabel = "x (μm)", ylabel = "y (μm)", zlabel = z_label,
               title  = "$(meta["sample"]) — $(meta["description"])\n$(sub_title)",
               titlegap = 12, titlesize = 14)
    sp = surface!(ax, xs, ys, h'; colormap = :viridis)
    Colorbar(fig[1, 2], sp, label = z_label)
    Label(fig[0, :], "drag = rotate · scroll = zoom · right-drag = pan · close window to exit";
          fontsize = 12, tellwidth = false)

    println("Loaded $(sample): $(M)×$(N), pixel_size = $(round(pixel_size, digits=4)) μm/px")
    println("h range: [$(round(minimum(h), digits=3)), $(round(maximum(h), digits=3))] " *
            (cal["method"] == "none" ? "(gauge-free)" : "μm"))

    screen = display(fig)
    println("\nInteractive window opened. Close it to exit.")
    wait(screen)
end

if abspath(PROGRAM_FILE) == @__FILE__
    sample = _parse_sample(ARGS)
    visualize(sample)
end
