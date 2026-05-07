# Reconstruct the Jun-series Vickers indent from `resources/Jun/pyramid/Q1..Q4.tif`.
# Outputs `data/jun_vickers_h.csv.gz` and `data/jun_vickers_meta.toml`.

include(joinpath(@__DIR__, "reconstruct_common.jl"))

const μmark_jun_vickers = 10.0 / (138 - 42)    # 0.1042 μm/px from SEM scale bar

const _DIR = joinpath(_PROJECT_ROOT, "resources", "Jun", "pyramid")

reconstruct_and_save(;
    sample          = "jun_vickers",
    description     = "Jun-series Vickers hardness indent (Jun/Q1..Q4) — square " *
                      "pyramidal impression made by a diamond indenter (136° face angle).",
    image_files     = [joinpath(_DIR, "Q$(i).tif") for i in 1:4],
    detector_angles = [0.0, π/2, π, 3π/2],
    pixel_size_um   = μmark_jun_vickers,
    calibration     = (method = :peak_to_trough, value = 6.0),
    # Vickers indent has very steep slopes near the apex — needs more
    # iterations than the 60-iteration default. `poly_order = 2` is on by
    # default (see scripts/reconstruct_common.jl) and is essential here.
    reconstruction_kwargs = (; max_iter = 200, tol = 1e-5),
)
