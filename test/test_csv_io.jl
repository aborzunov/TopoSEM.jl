@testset "csv_io — round-trip" begin
    include(joinpath(@__DIR__, "..", "scripts", "csv_io.jl"))

    Random.seed!(7)
    M, N = 64, 96
    h = randn(M, N) .* 12.5 .+ 0.7

    tmpdir = mktempdir()
    csv_path  = joinpath(tmpdir, "round_trip_h.csv.gz")
    toml_path = joinpath(tmpdir, "round_trip_meta.toml")

    write_topography_csv_gz(csv_path, h)
    h_back = read_topography_csv_gz(csv_path)

    @test size(h_back) == size(h)
    rel_err = maximum(abs.(h .- h_back)) / maximum(abs, h)
    @test rel_err < 1e-5

    meta = Dict(
        "sample" => "round_trip",
        "description" => "synthetic test fixture",
        "geometry" => Dict("pixel_size_um" => 0.0293,
                            "size_rows" => M,
                            "size_cols" => N),
        "calibration" => Dict("method" => "none"),
    )
    write_meta_toml(toml_path, meta)
    meta_back = read_meta_toml(toml_path)
    @test meta_back["sample"] == "round_trip"
    @test meta_back["geometry"]["pixel_size_um"] ≈ 0.0293
    @test meta_back["calibration"]["method"] == "none"

    rm(tmpdir; recursive = true, force = true)
end
