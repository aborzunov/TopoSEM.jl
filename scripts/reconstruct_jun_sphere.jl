# Reconstruct the Jun-series sphere/dome from `resources/Jun/sphere/S1..S4.tif`.
# Outputs `data/jun_sphere_h.csv.gz` and `data/jun_sphere_meta.toml`.

include(joinpath(@__DIR__, "reconstruct_common.jl"))

const μmark_jun_sphere = 20.0 / (85 - 41)      # 0.4545 μm/px from SEM scale bar

const _DIR = joinpath(_PROJECT_ROOT, "resources", "Jun", "sphere")

reconstruct_and_save(;
    sample          = "jun_sphere",
    description     = "Jun-series sphere / dome cap (Jun/S1..S4) — convex " *
                      "spherical cap of nominal 130 μm peak-to-trough height.",
    image_files     = [joinpath(_DIR, "S$(i).tif") for i in 1:4],
    detector_angles = [0.0, π/2, π, 3π/2],
    pixel_size_um   = μmark_jun_sphere,
    calibration     = (method = :peak_to_trough, value = 130.0),
    # IMPORTANT: explicitly override the package default `poly_order = 2`
    # back to 1 for the sphere. With poly_order = 2 the LSQ overfits the
    # smooth dome geometry and the reconstruction visually collapses into
    # a pyramid-like shape (see scripts/reconstruct_jun_sphere_exp2.jl
    # and the README "Polynomial-order observation" note).
    reconstruction_kwargs = (; poly_order = 1, max_iter = 200, tol = 1e-5),
)
