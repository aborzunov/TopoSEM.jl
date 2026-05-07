# Shared I/O helpers for the reproducibility pipeline.
#
# Topography matrices are persisted as gzipped CSV (rows × columns of Float64
# in `%.6e` format), and per-sample metadata as TOML sidecars. CSV is chosen
# for tooling agnosticism — pandas, NumPy, MATLAB, Excel, Julia all read it —
# and gzip wins back ~3× on disk relative to plain CSV.

using CodecZlib
using TOML
using Printf
using DelimitedFiles

const _CSV_FMT = Printf.Format("%.6e")

"""
    write_topography_csv_gz(path, h)

Write a 2-D `Float64` matrix `h` to a gzipped CSV at `path`. No header line;
one matrix row per text line; comma-separated; uses `%.6e` formatting
(≈ 7 significant digits, ≈ 13 chars per cell, ≈ 4 MB per 786 k cells after
gzip compression).
"""
function write_topography_csv_gz(path::AbstractString, h::AbstractMatrix)
    M, N = size(h)
    open(GzipCompressorStream, path, "w") do gz
        for i in 1:M
            for j in 1:N
                Printf.format(gz, _CSV_FMT, h[i, j])
                j < N && write(gz, ',')
            end
            write(gz, '\n')
        end
    end
    return path
end

"""
    read_topography_csv_gz(path) -> Matrix{Float64}

Load a topography produced by [`write_topography_csv_gz`](@ref).
"""
function read_topography_csv_gz(path::AbstractString)
    open(GzipDecompressorStream, path, "r") do gz
        return readdlm(gz, ',', Float64)
    end
end

"""
    write_meta_toml(path, meta::Dict)

Persist a metadata dictionary as TOML. The dict is expected to contain at
least the top-level keys `sample`, `description`, `source_files`, plus the
sub-tables `geometry`, `calibration`, `reconstruction`, `result`, `provenance`
(see `reconstruct_common.jl`).
"""
function write_meta_toml(path::AbstractString, meta::AbstractDict)
    open(path, "w") do io
        TOML.print(io, meta; sorted = false)
    end
    return path
end

"""
    read_meta_toml(path) -> Dict

Read back what [`write_meta_toml`](@ref) wrote. Sub-tables are accessed by
key (`meta["calibration"]["method"]`).
"""
function read_meta_toml(path::AbstractString)
    return TOML.parsefile(path)
end
