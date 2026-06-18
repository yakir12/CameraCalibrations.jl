using CameraCalibrations
using Test
using Aqua
using LinearAlgebra, Statistics
using FileIO, StaticArrays, PaddedViews, ColorTypes, LinearAlgebra, JSON, ImageTransformations

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
                      # piracies = (; treat_as_own = [Diagonal]),
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
            target = JSON.parsefile(joinpath(@__DIR__(), "example", "corners.json"))
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
            f = rectification(c, 1)
            i = RowCol(1,2)
            @test f(i) == c(i, 1)[[1,2]]
        end

    end

    @testset "Reference rectification" begin
        ei = 1
        sz = size(FileIO.load(CALIB.files[ei]))
        center = RowCol(sz[1] / 2, sz[2] / 2)
        north = RowCol(sz[1] / 4, sz[2] / 2)   # a quarter of the way down, i.e. above the center

        @testset "Invalid combinations error" begin
            @test_throws Exception only_scale(1.0; center = center)                    # center without index
            @test_throws Exception only_scale(1.0; north = north)                      # north without center/index
            @test_throws Exception only_scale(1.0; center = center, north = north)     # north + center without index
            @test_throws Exception fit(CALIB.files, CALIB.n_corners, CALIB.checker_size; center = center)
        end

        @testset "Valid combinations construct" begin
            @test only_scale(1.0) isa Calibration
            @test only_scale(1.0; extrinsic_index = 1) isa Calibration
            @test only_scale(1.0; extrinsic_index = 1, center = center) isa Calibration
            @test only_scale(1.0; extrinsic_index = 1, center = center, north = north) isa Calibration
        end

        c, _ = fit(CALIB.files, CALIB.n_corners, CALIB.checker_size; extrinsic_index = ei, center, north)

        @testset "Centering" begin
            @test rectification(c)(c.center) ≈ SVector(0, 0) atol = 1e-6
        end

        @testset "North alignment" begin
            p = rectification(c)(c.north)
            @test isapprox(p[1], 0; atol = 1e-6)
            @test p[2] > 0
        end

        @testset "Orientation (no mirror)" begin
            f = rectification(c)
            p0 = RowCol(sz[1] ÷ 2, sz[2] ÷ 2)
            ε = 1.0
            d1 = (f(p0 + RowCol(ε, 0)) - f(p0)) / ε
            d2 = (f(p0 + RowCol(0, ε)) - f(p0)) / ε
            J = hcat(d1, d2)
            @test det(J) > 0
        end

        @testset "Robustness" begin
            # extrinsic_index required
            cnone = fit(CALIB.files, CALIB.n_corners, CALIB.checker_size)[1]
            @test_throws Exception rectification(cnone)
            # center, no north → plain centered map (no rotation), and still centered
            ccenter, _ = fit(CALIB.files, CALIB.n_corners, CALIB.checker_size; extrinsic_index = ei, center)
            @test isnothing(ccenter.north)
            @test rectification(ccenter)(ccenter.center) ≈ SVector(0, 0) atol = 1e-6
            # index only → equivalent to the plain 2-arg deskew
            cidx, _ = fit(CALIB.files, CALIB.n_corners, CALIB.checker_size; extrinsic_index = ei)
            i = RowCol(1, 2)
            @test rectification(cidx)(i) == rectification(cidx, ei)(i)
            # existing 2-arg behaviour intact
            @test rectification(c, ei)(i) == c(i, ei)[[1, 2]]
            # out-of-range index on the direct-index entry points errors clearly (not a bare BoundsError)
            bad = length(c.extrinsics) + 1
            @test_throws ArgumentError rectification(c, bad)
            @test_throws ArgumentError c(i, bad)
            @test_throws ArgumentError c(CameraCalibrations.XYZ(0, 0, 0), bad)
            @test_throws ArgumentError rectification(c, 0)
        end

        @testset "JSON round-trip preserves metadata" begin
            copy = mktempdir() do path
                file = joinpath(path, "calibration.json")
                CameraCalibrations.save(file, c)
                CameraCalibrations.load(file)
            end
            @test copy.center ≈ c.center
            @test copy.north ≈ c.north
            @test copy.extrinsic_index == c.extrinsic_index
            @test rectification(copy)(copy.center) ≈ rectification(c)(c.center) atol = 1e-6
        end

        @testset "North too close to center" begin
            # below the 10-pixel minimum → error (angular direction is unreliable)
            @test_throws Exception fit(CALIB.files, CALIB.n_corners, CALIB.checker_size; extrinsic_index = ei, center, north = center .+ RowCol(5, 0))
            @test_throws Exception only_scale(1.0; extrinsic_index = 1, center = RowCol(100, 100), north = RowCol(103, 100))
            # exactly at the threshold (10 px) is allowed
            @test only_scale(1.0; extrinsic_index = 1, center = RowCol(100, 100), north = RowCol(110, 100)) isa Calibration
        end

        @testset "rectification factors through reference_affine" begin
            base = StaticArrays.pop ∘ c.image2real[c.extrinsic_index]
            f = CameraCalibrations.reference_affine(c) ∘ base
            i = RowCol(3, 4)
            @test rectification(c)(i) ≈ f(i)
        end

        @testset "Extrinsic index validation" begin
            # out-of-range index errors early (at construction), not later at use
            @test_throws Exception only_scale(1.0; extrinsic_index = 2, center = center)
            @test_throws Exception fit(CALIB.files, CALIB.n_corners, CALIB.checker_size; extrinsic_index = 0)
            @test_throws Exception fit(CALIB.files, CALIB.n_corners, CALIB.checker_size; extrinsic_index = length(CALIB.files) + 1)

            # select the extrinsic image by its file name (a String): exact tag, or by basename when tags are paths
            cstr, _ = fit(CALIB.files, CALIB.n_corners, CALIB.checker_size; extrinsic_index = CALIB.files[2])
            @test cstr.files[cstr.extrinsic_index] == CALIB.files[2]
            cbn, _ = fit(CALIB.files, CALIB.n_corners, CALIB.checker_size; extrinsic_index = basename(CALIB.files[2]))
            @test cbn.files[cbn.extrinsic_index] == CALIB.files[2]
            @test_throws Exception fit(CALIB.files, CALIB.n_corners, CALIB.checker_size; extrinsic_index = "not-a-file.png")
            # an ambiguous basename (two images share it) errors rather than picking one
            dup_tags = ["dirA/x.png", "dirB/x.png"]
            dup_imgs = [Gray.(FileIO.load(CALIB.files[1])), Gray.(FileIO.load(CALIB.files[2]))]
            @test_throws Exception fit(dup_tags, dup_imgs, CALIB.n_corners, CALIB.checker_size; extrinsic_index = "x.png")

            # a dropped image (failed corner detection) must not shift which image extrinsic_index refers to
            imgs = [Gray.(FileIO.load(f)) for f in CALIB.files]
            blank = fill(zero(eltype(imgs[1])), size(imgs[1]))   # detection fails here → image dropped
            tags2 = ["blank", basename.(CALIB.files)...]
            imgs2 = [blank, imgs...]
            target = basename(CALIB.files[2])                    # sits after the blank in the original list
            pos = findfirst(==(target), tags2)
            cdrop, _ = fit(tags2, imgs2, CALIB.n_corners, CALIB.checker_size; extrinsic_index = pos)
            @test !("blank" in cdrop.files)                      # the blank really was dropped
            @test cdrop.files[cdrop.extrinsic_index] == target   # index still points to the intended image
            # choosing the undetectable image is an error (the calibration is unusable as a reference)
            @test_throws Exception fit(tags2, imgs2, CALIB.n_corners, CALIB.checker_size; extrinsic_index = "blank")
        end
    end

    @testset "Lens distortion inverse" begin
        # round-trip forward → inverse must recover the undistorted point, including barrel (negative k1) lenses,
        # over the well-posed radius range. (The previous Optim-based inverse was inaccurate for barrel distortion.)
        ld    = CameraCalibrations.lens_distortion
        inv_ld = CameraCalibrations.inv_lens_distortion
        for k in ([-0.2, 0.0, 0.0], [-0.3, 0.1, 0.0], [-0.3, 0.1, -0.05], [0.2, -0.05, 0.0], [0.4, 0.0, 0.0])
            for r in 0.0:0.1:0.9, θ in range(0, 2π, length = 9)[1:8]
                v  = SVector(r*cos(θ), r*sin(θ))
                @test inv_ld(ld(v, k), k) ≈ v atol = 1e-8
            end
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
                # `copy` is always a `Calibration` (load reconstructs one), so its metadata is
                # narrowed to `Float32`; compare with a tolerance loose enough for that narrowing.
                @test isapprox(org.center, copy.center; atol = 1e-5)
                @test isapprox(org.north, copy.north; atol = 1e-5)
                @test org.extrinsic_index == copy.extrinsic_index
            end
        end
        @testset "no metadata round-trips" begin
            c0 = only_scale(1.0)   # center / north / extrinsic_index all nothing
            copy0 = mktempdir() do path
                file = joinpath(path, "c.json")
                CameraCalibrations.save(file, c0)
                CameraCalibrations.load(file)
            end
            @test isnothing(copy0.center)
            @test isnothing(copy0.north)
            @test isnothing(copy0.extrinsic_index)
        end

        @testset "invalid CalibrationIO cannot be built" begin
            # validation happens up front, so an unsaveable/unloadable combination errors immediately
            mk(; kw...) = CameraCalibrations.CalibrationIO(rand(2, 2), [rand(3, 3)], rand(3, 3), rand(3), ["a"]; kw...)
            @test_throws Exception mk(center = [1.0, 2.0])                      # center without index
            @test_throws Exception mk(extrinsic_index = 1, center = [1.0, 2.0], north = [1.5, 2.0])  # north < 10 px
        end

        @testset "MATLAB" begin
            matfile = joinpath(@__DIR__(), "example", "cameraParams.mat")
            @test_nowarn CameraCalibrations.loadMAT(matfile)
            # explicit loaders reject the wrong format with a helpful error
            @test_throws Exception CameraCalibrations.load(matfile)            # MATLAB file through the JSON loader
            mktempdir() do path
                file = joinpath(path, "c.json")
                CameraCalibrations.save(file, only_scale(1.0))
                @test_throws Exception CameraCalibrations.loadMAT(file)        # JSON file through the MATLAB loader
            end
        end

        @testset "only scale" begin
            for factor in (0.5, 1, 2)
                c = only_scale(factor)
                calib = rectification(c, 1)
                rc = RowCol(1, 2)
                @test all(calib(rc) .≈ rc/factor)
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

        @testset "Rectified extrinsic" begin
            mktempdir() do plot_folder
                sz = size(FileIO.load(CALIB.files[1]))
                center = RowCol(sz[1] / 2, sz[2] / 2)
                north = RowCol(sz[1] / 4, sz[2] / 2)
                fit(CALIB.files, CALIB.n_corners, CALIB.checker_size; CALIB.aspect, plot_folder, extrinsic_index = 1, center, north)
                @test "rectified_extrinsic.png" in readdir(plot_folder)
            end
        end

        @testset "In-memory images plot from memory (no reload from disk)" begin
            # tags are not file paths, so plotting must use the in-memory images, not FileIO.load(tag)
            imgs = [Gray.(FileIO.load(f)) for f in CALIB.files]
            tags = ["frame$i" for i in eachindex(CALIB.files)]   # bare tags, no extension
            mktempdir() do plot_folder
                sz = size(imgs[1])
                center = RowCol(sz[1] / 2, sz[2] / 2)
                north = RowCol(sz[1] / 4, sz[2] / 2)
                fit(tags, imgs, CALIB.n_corners, CALIB.checker_size; plot_folder, extrinsic_index = 1, center, north)
                @test "rectified_extrinsic.png" in readdir(plot_folder)
                @test "frame1.png" in readdir(plot_folder)        # calibration plot, extension appended
            end
        end
    end

end
