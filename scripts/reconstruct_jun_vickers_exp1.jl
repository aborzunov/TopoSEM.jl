# Experiment #1 for the Jun-series Vickers indent — same input data as
# jun_vickers but reconstructed with `poly_order = 1` (linear φ).
#
# Preserved as a CAUTIONARY example: the linear model cannot represent
# the indent's near-vertical wall slopes, the iteration plateaus at
# `‖δh‖ ≈ 100` and `‖δf‖ ≈ 15000` regardless of `max_iter`, and the wall
# detail is visibly under-fit. The companion `jun_vickers` script keeps
# `poly_order = 2`, which captures the walls cleanly.

include(joinpath(@__DIR__, "reconstruct_common.jl"))

const μmark_jun_vickers = 10.0 / (138 - 42)    # 0.1042 μm/px

const _DIR = joinpath(_PROJECT_ROOT, "resources", "Jun", "pyramid")

reconstruct_and_save(;
    sample          = "jun_vickers_exp1",
    description     = "Jun-series Vickers indent — EXPERIMENT №1: reconstructed " *
                      "with poly_order = 1 (linear φ). The indent walls are " *
                      "essentially vertical and the linear model cannot " *
                      "represent them; preserved as a contrast to the working " *
                      "poly_order = 2 result in jun_vickers.",
    image_files     = [joinpath(_DIR, "Q$(i).tif") for i in 1:4],
    detector_angles = [0.0, π/2, π, 3π/2],
    pixel_size_um   = μmark_jun_vickers,
    calibration     = (method = :peak_to_trough, value = 6.0),
    reconstruction_kwargs = (; poly_order = 1, max_iter = 200, tol = 1e-5),
)
