# Reconstruct the May calibration sphere from `resources/May/{13,16,1,4}.jpg`.
# This is the May-series sphere acquisition that SurfaceTopography.jl uses
# as a 2D-apparatus-function calibration target. Channel order (A,B,C,D)
# is identified by the SurfaceTopography sample list (see survey notes).

include(joinpath(@__DIR__, "reconstruct_common.jl"))

const μmark_may_sphere = 100.0 / (1006 - 698)  # 0.3247 μm/px from SEM scale bar

const _DIR = joinpath(_PROJECT_ROOT, "resources", "May")

reconstruct_and_save(;
    sample          = "may_sphere",
    description     = "May calibration sphere — 4 BSE quadrants " *
                      "(A=13.jpg, B=16.jpg, C=1.jpg, D=4.jpg). " *
                      "Used as the apparatus-function calibration target in " *
                      "the companion SurfaceTopography.jl pipeline.",
    image_files     = [joinpath(_DIR, "$(idx).jpg") for idx in (13, 16, 1, 4)],
    detector_angles = [0.0, π/2, π, 3π/2],
    pixel_size_um   = μmark_may_sphere,
    calibration     = (method = :peak_to_trough, value = 100.0),
)
