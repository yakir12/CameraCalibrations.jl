using CameraCalibrations
using Test
using Aqua
using LinearAlgebra, Statistics
using FileIO, StaticArrays, PaddedViews, ColorTypes, LinearAlgebra, JSON3, ImageTransformations

files = filter(endswith(".png"), readdir(joinpath(@__DIR__(), "example"), join = true))
n_corners = (5, 8)
checker_size = 1
aspect = 1
# plot_folder
const CALIB = (; files, n_corners, checker_size, aspect)

function index2bw(ij::CartesianIndex) 
    i, j = Tuple(ij)
    isodd(i) ? !isodd(j) : isodd(j)
end

function generate_checkerboard(n_corners, n) 
    xys = [n .* SVector{2, Float32}(Tuple(ij)) - SVector(0.5, 0.5) for ij in CartesianIndices(StepRange.(2, 1, n_corners .+ 1))]
    # reverse!(xys, dims = issorted(n_corners) ? 1 : 2)
    img = index2bw.(CartesianIndices(n_corners .+ 1))
    imgl = kron(PaddedView(true, img, UnitRange.(0, n_corners .+ 2)), ones(Int, n, n))
    return xys, Gray.(imgl)
end

function calc_rms(xys1::Matrix, xys2::Matrix)
    min(sqrt(mean(LinearAlgebra.norm_sqr.(xys2 .- xys1))), sqrt(mean(LinearAlgebra.norm_sqr.(reverse(xys2) .- xys1))))
end

function calc_rms(n_corners::NTuple{2, Int}, ratio::Int)
    xys, img = generate_checkerboard(n_corners, ratio)
    res = CameraCalibrations._detect_corners(img, n_corners)
    return !ismissing(res) && calc_rms(xys, res)
end

@testset "CameraCalibrations.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(CameraCalibrations; 
                      piracies = (; treat_as_own = [Diagonal]),
                      stale_deps = (; ignore = [:ImageIO]))
    end

    @testset "Detect corners" begin

        @testset "In artificial images" begin
            for w in 13:15, h in 13:15, ratio in 95:100
                if isodd(w) ≠ isodd(h)
                    n_corners = (w, h)
                    @test calc_rms(n_corners, ratio) ≈ 0
                end
            end
        end

        if Threads.nthreads() > 1
            @testset "In artificial images (threaded)" begin
                xs = [((w, h), r) for w in 13:15 for h in 13:15 for r in 95:100 if isodd(w) ≠ isodd(h)]
                Threads.@threads for i in 1:length(xs)
                    n_corners, ratio = xs[i]
                    @test calc_rms(n_corners, ratio) ≈ 0
                end
            end
        else
            @info "Julia is running with only 1 thread, skipping multithreaded tests"
        end

        @testset "In real images" begin
            target = JSON3.read(joinpath(@__DIR__(), "example", "corners.json"))
            for file in CALIB.files
                img = FileIO.load(file)
                res = CameraCalibrations._detect_corners(img, CALIB.n_corners)
                target_corners = target[basename(file)]
                @test !ismissing(res) && calc_rms(collect(reshape(target_corners, CALIB.n_corners)), res) < 1e-4
            end
        end

    end

    @testset "Full calibration" begin

        c, (n, ϵ...) = fit(CALIB.files, CALIB.n_corners, CALIB.checker_size)

        @testset "$k accuracy" for (k, v) in pairs(ϵ)
            @test v < 1
        end

        @testset "Rectification" begin
            extrinsic_index = 1
            extrinsic_file = files[extrinsic_index]
            f = rectification(c, 1)
            i = RowCol(1,2)
            @test f(i) == c(i, 1)[[1,2]]
        end

    end

    @testset "IO" begin
        @testset "$c" for c in (Calibration, CalibrationIO)
            for _ in 1:100
                org = rand(c)
                copy = mktempdir() do path
                    file = joinpath(path, "calibration.json")
                    CameraCalibrations.save(file, org)
                    CameraCalibrations.load(file)
                end
                @test org.files == copy.files
            end
        end
    end

    @testset "Full calibration with a different $aspect ratio" for aspect in (0.75, 1.25)
        mktempdir() do path
            for file in CALIB.files
                img = FileIO.load(file)
                img = imresize(img; ratio = (1, aspect))
                FileIO.save(joinpath(path, basename(file)), img)
            end

            files = readdir(path, join = true)
            c, (n, ϵ...) = fit(files, CALIB.n_corners, CALIB.checker_size; aspect)

            @testset "$k accuracy" for (k, v) in pairs(ϵ)
                @test v < 1
            end

            @testset "Rectification" begin
                extrinsic_index = 1
                extrinsic_file = files[extrinsic_index]
                f = rectification(c, 1)
                i = RowCol(1,2)
                @test f(i) == c(i, 1)[[1,2]]
            end
        end
    end

    if Threads.nthreads() > 1
        @testset "Full calibration (threaded)" begin 
            Threads.@threads for _ in 1:5
                c, (n, ϵ...) = fit(CALIB.files, CALIB.n_corners, CALIB.checker_size; CALIB.aspect)
                @testset "$k accuracy" for (k, v) in pairs(ϵ)
                    @test v < 1
                end
            end
        end
    else
        @info "Julia is running with only 1 thread, skipping multithreaded tests"
    end

    @testset "Plotting" begin 
        mktempdir() do plot_folder
            fit(CALIB.files, CALIB.n_corners, CALIB.checker_size; CALIB.aspect, plot_folder)
            @test readdir(plot_folder) == basename.(CALIB.files)
        end
    end

end
