# Reconstruct the May piezoceramic sample (microcracked surface, "bubbles")
# from `resources/May/{34,35,26,27}.jpg` (channels A, B, C, D).

include(joinpath(@__DIR__, "reconstruct_common.jl"))

const μmark_may_piezoceramic = 10.0 / (939 - 768)   # 0.0585 μm/px

const _DIR = joinpath(_PROJECT_ROOT, "resources", "May")

reconstruct_and_save(;
    sample          = "may_piezoceramic",
    description     = "May piezoceramic sample — microcracked surface, " *
                      "4 BSE quadrants (A=34.jpg, B=35.jpg, C=26.jpg, D=27.jpg).",
    image_files     = [joinpath(_DIR, "$(idx).jpg") for idx in (34, 35, 26, 27)],
    detector_angles = [0.0, π/2, π, 3π/2],
    pixel_size_um   = μmark_may_piezoceramic,
    calibration     = nothing,                       # gauge-free
)
