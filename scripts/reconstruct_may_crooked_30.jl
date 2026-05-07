# Reconstruct the May "crooked pyramid" at 30° tilt (companion to the 60°
# variant). Same sample, different SEM detector tilt — useful for
# comparing model fidelity across acquisition geometries.

include(joinpath(@__DIR__, "reconstruct_common.jl"))

const μmark_may_crooked_pyramid = 30.0 / (995 - 687)   # 0.0974 μm/px

const _DIR = joinpath(_PROJECT_ROOT, "resources", "May")

reconstruct_and_save(;
    sample          = "may_crooked_30",
    description     = "May crooked-pyramid sample at 30° detector tilt, " *
                      "4 BSE quadrants (A=46.jpg, B=47.jpg, C=40.jpg, D=41.jpg).",
    image_files     = [joinpath(_DIR, "$(idx).jpg") for idx in (46, 47, 40, 41)],
    detector_angles = [0.0, π/2, π, 3π/2],
    pixel_size_um   = μmark_may_crooked_pyramid,
    calibration     = nothing,
)
