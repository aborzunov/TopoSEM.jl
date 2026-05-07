# Reconstruct the May "crooked pyramid" (tilted Vickers indent) at 60° tilt
# from `resources/May/{45,48,39,42}.jpg` (channels A, B, C, D).

include(joinpath(@__DIR__, "reconstruct_common.jl"))

const μmark_may_crooked_pyramid = 30.0 / (995 - 687)   # 0.0974 μm/px

const _DIR = joinpath(_PROJECT_ROOT, "resources", "May")

reconstruct_and_save(;
    sample          = "may_crooked_60",
    description     = "May crooked-pyramid sample at 60° detector tilt, " *
                      "4 BSE quadrants (A=45.jpg, B=48.jpg, C=39.jpg, D=42.jpg).",
    image_files     = [joinpath(_DIR, "$(idx).jpg") for idx in (45, 48, 39, 42)],
    detector_angles = [0.0, π/2, π, 3π/2],
    pixel_size_um   = μmark_may_crooked_pyramid,
    calibration     = nothing,
)
