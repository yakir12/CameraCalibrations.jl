# using Revise
using CameraCalibrations
using Test
using Aqua
using FileIO, StaticArrays, PaddedViews, ColorTypes, LinearAlgebra, JSON3, ImageTransformations

function index2bw(ij::CartesianIndex) 
    i, j = Tuple(ij)
    isodd(i) ? isodd(j) : !isodd(j)
end

function generate_checkerboard(n_corners, n) 
    xys = [n .* SVector{2, Float32}(Tuple(ij)) - SVector(0.5, 0.5) for ij in CartesianIndices(StepRange.(n_corners .+ 1, -1, 2))]
    # reverse!(xys, dims = 1)
    img = index2bw.(CartesianIndices(n_corners .+ 1))
    imgl = kron(PaddedView(true, img, UnitRange.(0, n_corners .+ 2)), ones(Int, n, n))
    # imgw = imfilter(imgl, Kernel.gaussian(2))
    return xys, imgl
end

function compare(n_corners, ratio)
    xys, img = generate_checkerboard(n_corners, ratio)
    mktempdir() do path
        # path = mktempdir(; cleanup = false)
        file = joinpath(path, "img.png")
        FileIO.save(file, Gray.(img))
        res = CameraCalibrations._detect_corners(file, n_corners)
        return !ismissing(res) && all(<(1), norm.(xys .- res[2]))
    end
end

# using ImageDraw
# dir = joinpath(@__DIR__(), "example")
# files = filter(file -> last(splitext(file)) == ".png", readdir(dir, join = true))
# n = length(files)
# n_corners = (5, 8)
# d = Dict()
# for file in files
#     res = CameraCalibrations._detect_corners(file, n_corners)
#     if !ismissing(res) 
#         _, xys = res
#         img = RGB.(FileIO.load(file))
#         for xy in xys
#             draw!(img, Cross(Point(round.(Int, reverse(xy))...), 5), RGB(1,0,0))
#         end
#         FileIO.save(basename(file), img)
#         d[basename(file)] = vec(xys)
#     end
# end
# open("example/corners.json", "w") do io
#     print(io, JSON3.write(d))
# end

@testset "CameraCalibrations.jl" begin
    @testset "Code quality (Aqua.jl)" begin
        Aqua.test_all(CameraCalibrations; piracies = false)
    end

    @testset "Detect corners" begin
        @testset "In artificial images" begin
            for w in 13:15, h in 13:15, ratio in 95:100
                if isodd(w) ≠ isodd(h)
                    w, h, ratio = 14, 13, 100
                    n_corners = (w, h)
                    @test compare(n_corners, ratio)
                end
            end
        end

        # # Threading doesn't work at all
        # @testset "In artificial images (threaded)" begin
        #     xs = [((w, h), r) for w in 13:15 for h in 13:15 for r in 95:100 if isodd(w) ≠ isodd(h)]
        #     Threads.@threads for i in 1:length(xs)
        #         n_corners, ratio = xs[i]
        #         @test compare(n_corners, ratio)
        #     end
        # end

        @testset "In real images" begin
            dir = joinpath(@__DIR__(), "example")
            files = filter(file -> last(splitext(file)) == ".png", readdir(dir, join = true))
            n = length(files)
            n_corners = (5, 8)
            target = JSON3.read(joinpath(dir, "corners.json"))
            for file in files
                res = CameraCalibrations._detect_corners(file, n_corners)
                target_corners = target[basename(file)]
                @test !ismissing(res) && all(<(1), norm.(target_corners .- vec(res[2])))
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
                img = imresize(img; ratio = (aspect, 1))
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


    end
