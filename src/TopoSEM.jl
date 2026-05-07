module TopoSEM

using LinearAlgebra
using Statistics
using Printf
using StaticArrays
using OffsetArrays
using FFTW
using TiffImages
using FileIO
using ImageIO       # registers PNG/etc. backends (transitive)
using JpegTurbo     # registers JPEG backend with FileIO
using ColorTypes: Gray

include("poly.jl")
include("detector.jl")
include("derivatives.jl")
include("weights.jl")
include("step_a.jl")
include("step_b.jl")
include("step_c.jl")
include("main.jl")
include("io.jl")

export BivariatePoly, DetectorModel, TopoSEMResult,
       forward_model, forward_image,
       gradient_centered, laplacian_5pt,
       pixel_weights,
       solve_pixelwise!, integrate_fourier, recalibrate!,
       topo_sem,
       load_bse_image, load_bse_images, load_bse_quadrants, detect_data_rows,
       detrend_planar, calibrate_height

end # module
