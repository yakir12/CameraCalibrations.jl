using CameraCalibrations
using Test
using Aqua
using LinearAlgebra, Statistics
using FileIO, StaticArrays, PaddedViews, ColorTypes, LinearAlgebra, JSON3, ImageTransformations

function index2bw(ij::CartesianIndex) 
    i, j = Tuple(ij)
    isodd(i) ? !isodd(j) : isodd(j)
end

function generate_checkerboard(n_corners, n) 
    xys = [n .* SVector{2, Float32}(Tuple(ij)) - SVector(0.5, 0.5) for ij in CartesianIndices(StepRange.(2, 1, n_corners .+ 1))]
    # reverse!(xys, dims = issorted(n_corners) ? 1 : 2)
    img = index2bw.(CartesianIndices(n_corners .+ 1))
    imgl = kron(PaddedView(true, img, UnitRange.(0, n_corners .+ 2)), ones(Int, n, n))
    return xys, imgl
end

function calc_rms(xys1::Matrix, xys2::Matrix)
    min(sqrt(mean(LinearAlgebra.norm_sqr.(xys2 .- xys1))), sqrt(mean(LinearAlgebra.norm_sqr.(reverse(xys2) .- xys1))))
end

function calc_rms(n_corners::NTuple{2, Int}, ratio::Int)
    xys, img = generate_checkerboard(n_corners, ratio)
    mktempdir() do path
        file = joinpath(path, "img.png")
        FileIO.save(file, Gray.(img))
        res = CameraCalibrations._detect_corners(file, n_corners)
        return !ismissing(res) && calc_rms(xys, res[2])
    end
end

@testset "CameraCalibrations.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(CameraCalibrations; piracies = false)
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

        @testset "In artificial images (threaded)" begin
            xs = [((w, h), r) for w in 13:15 for h in 13:15 for r in 95:100 if isodd(w) ≠ isodd(h)]
            Threads.@threads for i in 1:length(xs)
                n_corners, ratio = xs[i]
                @test calc_rms(n_corners, ratio) ≈ 0
            end
        end

        @testset "In real images" begin
            dir = joinpath(@__DIR__(), "example")
            files = filter(file -> last(splitext(file)) == ".png", readdir(dir, join = true))
            n = length(files)
            n_corners = (5, 8)
            target = JSON3.read(joinpath(dir, "corners.json"))
            for file in files
                res = CameraCalibrations._detect_corners(file, n_corners)
                target_corners = target[basename(file)]
                @test !ismissing(res) && calc_rms(collect(reshape(target_corners, n_corners)), res[2]) < 1e-4
            end
        end

    end

    @testset "Full calibration" begin

        n_corners = (5, 8)
        dir = joinpath(@__DIR__(), "example")
        files = filter(file -> last(splitext(file)) == ".png", readdir(dir, join = true))
        checker_size = 1
        c, (n, ϵ...) = fit(files, n_corners, checker_size)

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
        for _ in 1:100
            org = rand(CameraCalibrations.CalibrationIO)
            copy = mktempdir() do path
                file = joinpath(path, "calibration.json")
                CameraCalibrations.save(file, org)
                CameraCalibrations.load(file)
            end
            @test org.files == copy.files
        end
    end

    @testset "Full calibration with a different $aspect ratio" for aspect in (0.75, 1.25)
        n_corners = (5, 8)
        dir = joinpath(@__DIR__(), "example")
        files = filter(file -> last(splitext(file)) == ".png", readdir(dir, join = true))
        checker_size = 1
        mktempdir() do path
            for file in files
                img = FileIO.load(file)
                img = imresize(img; ratio = (1, aspect))
                FileIO.save(joinpath(path, basename(file)), img)
            end

            files = readdir(path, join = true)
            c, (n, ϵ...) = fit(files, n_corners, checker_size; aspect)

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

    @testset "Full calibration with different aspect ratios (threaded)" begin 
        aspects = 0.75:0.1:1.25
        n_corners = (5, 8)
        dir = joinpath(@__DIR__(), "example")
        files = filter(file -> last(splitext(file)) == ".png", readdir(dir, join = true))
        checker_size = 1
        Threads.@threads for i in 1:length(aspects)
            mktempdir() do path
                for file in files
                    img = FileIO.load(file)
                    img = imresize(img; ratio = (1, aspects[i]))
                    FileIO.save(joinpath(path, basename(file)), img)
                end

                files = readdir(path, join = true)
                c, (n, ϵ...) = fit(files, n_corners, checker_size; aspect = aspects[i])

                @testset "$k accuracy" for (k, v) in pairs(ϵ)
                    @test v < 1
                end

            end
        end
    end

end
