# Reconstruct Feb sample group P (Feb/PA, PB, PC, PD) — small pyramid feature.

include(joinpath(@__DIR__, "reconstruct_common.jl"))

const μmark_feb_pyramid = 10.0 / (955 - 752)   # 0.0493 μm/px from SEM scale bar

const _DIR = joinpath(_PROJECT_ROOT, "resources", "Feb")

reconstruct_and_save(;
    sample          = "feb_P",
    description     = "Feb pyramid sample (PA..PD), 4-quadrant BSE.",
    image_files     = [joinpath(_DIR, "P$(c).tif") for c in "ABCD"],
    detector_angles = [0.0, π/2, π, 3π/2],
    pixel_size_um   = μmark_feb_pyramid,
    calibration     = nothing,                  # gauge-free
)
