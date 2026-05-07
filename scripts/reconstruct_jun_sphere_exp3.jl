# Experiment #3 for the Jun-series sphere — same input data, but with
# `poly_order = 3` (cubic Taylor expansion of φ).
#
# README (§ "Polynomial-order observation") notes that orders ≥ 3 are
# expected to be unstable without ridge regularisation in step C — this
# script preserves the actual outcome on the dome for direct comparison.

include(joinpath(@__DIR__, "reconstruct_common.jl"))

const μmark_jun_sphere = 20.0 / (85 - 41)      # 0.4545 μm/px

const _DIR = joinpath(_PROJECT_ROOT, "resources", "Jun", "sphere")

reconstruct_and_save(;
    sample          = "jun_sphere_exp3",
    description     = "Jun-series sphere — EXPERIMENT №3 (poly_order = 3, " *
                      "cubic φ). Probes whether higher polynomial degree " *
                      "improves dome reconstruction over the working " *
                      "poly_order = 1 baseline.",
    image_files     = [joinpath(_DIR, "S$(i).tif") for i in 1:4],
    detector_angles = [0.0, π/2, π, 3π/2],
    pixel_size_um   = μmark_jun_sphere,
    calibration     = (method = :peak_to_trough, value = 130.0),
    reconstruction_kwargs = (; poly_order = 3, max_iter = 200, tol = 1e-5),
)
