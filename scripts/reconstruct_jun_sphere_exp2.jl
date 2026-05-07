# Experiment #2 for the Jun-series sphere — same input data as jun_sphere
# but reconstructed with `poly_order = 2`.
#
# This is preserved as a CAUTIONARY example: on the dome the quadratic φ
# overfits and the reconstruction visually collapses into a pyramid-like
# shape with asymmetric h-range (-38.6 .. +91.4 μm). The companion
# `jun_sphere` script keeps `poly_order = 1`, which gives the correct
# convex-cap reconstruction.

include(joinpath(@__DIR__, "reconstruct_common.jl"))

const μmark_jun_sphere = 20.0 / (85 - 41)      # 0.4545 μm/px

const _DIR = joinpath(_PROJECT_ROOT, "resources", "Jun", "sphere")

reconstruct_and_save(;
    sample          = "jun_sphere_exp2",
    description     = "Jun-series sphere — EXPERIMENT №2: reconstructed with " *
                      "poly_order = 2 (the new package default). On this " *
                      "dome the quadratic φ over-fits and the reconstruction " *
                      "visually resembles a pyramid; preserved as a contrast " *
                      "to the working poly_order = 1 result in jun_sphere.",
    image_files     = [joinpath(_DIR, "S$(i).tif") for i in 1:4],
    detector_angles = [0.0, π/2, π, 3π/2],
    pixel_size_um   = μmark_jun_sphere,
    calibration     = (method = :peak_to_trough, value = 130.0),
    reconstruction_kwargs = (; poly_order = 2, max_iter = 200, tol = 1e-5),
)
