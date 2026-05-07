# Reconstruct Feb sample group S (Feb/SA..SD) — sphere feature.

include(joinpath(@__DIR__, "reconstruct_common.jl"))

const μmark_feb_sphere = 100.0 / (990 - 715)   # 0.3636 μm/px

const _DIR = joinpath(_PROJECT_ROOT, "resources", "Feb")

reconstruct_and_save(;
    sample          = "feb_S",
    description     = "Feb sphere sample (SA..SD), 4-quadrant BSE.",
    image_files     = [joinpath(_DIR, "S$(c).tif") for c in "ABCD"],
    detector_angles = [0.0, π/2, π, 3π/2],
    pixel_size_um   = μmark_feb_sphere,
    calibration     = nothing,
)
