# Reconstruct Feb sample group PR (Feb/PRA..PRD) — "crooked pyramid".

include(joinpath(@__DIR__, "reconstruct_common.jl"))

const μmark_crooked_pyramid = 30.0 / (995 - 687)   # 0.0974 μm/px

const _DIR = joinpath(_PROJECT_ROOT, "resources", "Feb")

reconstruct_and_save(;
    sample          = "feb_PR",
    description     = "Feb 'crooked pyramid' sample (PRA..PRD), rotated/skewed " *
                      "pyramidal feature, 4-quadrant BSE.",
    image_files     = [joinpath(_DIR, "PR$(c).tif") for c in "ABCD"],
    detector_angles = [0.0, π/2, π, 3π/2],
    pixel_size_um   = μmark_crooked_pyramid,
    calibration     = nothing,
)
