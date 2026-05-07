# Reconstruct Feb sample group V (Feb/VA..VD) — V-groove feature.
# Same SEM acquisition as Feb/P, so same pixel scale.

include(joinpath(@__DIR__, "reconstruct_common.jl"))

const μmark_feb_V = 10.0 / (955 - 752)         # 0.0493 μm/px (same as feb_P)

const _DIR = joinpath(_PROJECT_ROOT, "resources", "Feb")

reconstruct_and_save(;
    sample          = "feb_V",
    description     = "Feb V-groove sample (VA..VD), 4-quadrant BSE.",
    image_files     = [joinpath(_DIR, "V$(c).tif") for c in "ABCD"],
    detector_angles = [0.0, π/2, π, 3π/2],
    pixel_size_um   = μmark_feb_V,
    calibration     = nothing,
)
