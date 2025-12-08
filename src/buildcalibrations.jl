function obj2img(Rs, ts, frow, fcol, crow, ccol, checker_size)
    intrinsic = AffineMap(SDiagonal(frow, fcol), SVector(crow, ccol))
    extrinsics = AffineMap.(Base.splat(RotationVec).(Rs), SVector{3, Float64}.(ts))
    scale = LinearMap(SDiagonal{3}(I/checker_size))
    return intrinsic, extrinsics, scale
end

"""
    Calibration(tags::Vector{AbstractString}, imgs::Vector{Matrix{Gray}}, n_corners, checker_size; aspect = 1, radial_parameters = 1, plot_folder::Union{Nothing, String} = nothing)

Build a calibration object. `tags` are references (just names) to each of the images in `imgs`. Typically, one of these will be the chosen image for the extrinsics parameters. `imgs` are the images of the checkerboard. `n_corners` is a tuple of the number of corners in each of the sides of the checkerboard. `checker_size` is the physical size of the checker (e.g. in cm). `aspect` is the aspect ratio of the images in `imgs` (most commonly 1). `radial_parameters` controls how many radial lens distortion parameters are included in the model. Finally, `plot_folder` is the path to a directory where diagnostic images can be saved to.
"""
function fit(tags::Vector{T}, imgs::Vector{Matrix{S}}, n_corners, checker_size; aspect = 1, radial_parameters::Int = 1, plot_folder::Union{Nothing, String} = nothing) where {T <: AbstractString, S <: Gray}
    @assert length(tags) == length(imgs) "`tags` and `imgs` should have the same length"
    files, objpoints, imgpointss, sz, k, Rs, ts, frow, fcol, crow, ccol = detect_fit(tags, imgs, n_corners, radial_parameters, aspect)
    objpoints = objpoints .* checker_size
    intrinsic, extrinsics, scale = obj2img(Rs, ts, frow, fcol, crow, ccol, checker_size)
    c = Calibration(intrinsic, extrinsics, scale, k, files)

    # cf = improve(c, improve_n, improve_threshold, aspect, with_distortion)

    ϵ = calculate_errors(c, imgpointss, objpoints, checker_size, sz, files, n_corners)
    plot(plot_folder, c, imgpointss, n_corners, checker_size, sz)
    return (c, ϵ)
end

"""
    fit(files::Vector{AbstractString}, n_corners, checker_size; aspect = 1, radial_parameters = 1, plot_folder::Union{Nothing, String} = nothing)

Build a calibration object. `files` are file names to the images of the checkerboard.
"""
function fit(files::Vector{T}, n_corners, checker_size; kwargs...) where T <: AbstractString
    imgs = [Gray.(FileIO.load(file)) for file in files]
    fit(files, imgs, n_corners, checker_size; kwargs...)
end

function _reprojection(c, i, objpoints, imgpoints)
    reprojected = c.(objpoints, i)
    sum(LinearAlgebra.norm_sqr, reprojected .- imgpoints)
end

"""
    calculate_errors(c)
Calculate reprojection, projection, distance, and inverse errors for the calibration `c`. `distance` measures the mean error of the distance between all adjacent checkerboard corners from the expected `checker_size`. `inverse` measures the mean error of applying the calibration's transformation and its inverse `inverse_samples` times.
"""
function calculate_errors(c, imgpointss, objpoints, checker_size, sz, files, n_corners)
    inverse_samples = 100
    reprojection = 0.0
    projection = 0.0
    distance = 0.0
    inverse = 0.0
    for (i, imgpoints) in pairs(imgpointss)
        reprojection += _reprojection(c, i, objpoints, imgpoints)

        projected = c.(imgpoints, i)
        projection += sum(LinearAlgebra.norm_sqr, projected .- objpoints)

        distance += sum(1:2) do dims
            sum(abs2, norm.(diff(projected; dims)) .- checker_size)
        end

        inverse += sum(1:inverse_samples) do _
            rc = rand(RowCol) .* (sz .- 1) .+ 1
            projected = c(rc, i)
            reprojected = c(projected, i)
            LinearAlgebra.norm_sqr(rc .- reprojected)
        end
    end
    n_files = length(files)
    n = prod(n_corners)*n_files
    reprojection = sqrt(reprojection/n)
    projection = sqrt(projection/n)
    distance = sqrt(distance/prod(n_corners .- 1)/n_files)
    inverse = sqrt(inverse/inverse_samples/n_files)
    return (; n = n_files, reprojection, projection, distance, inverse)
end

function only_scale(scale_factor)
    intrinsic = AffineMap{Matrix{Float64}, Vector{Float64}}(collect(I(2)), zeros(2))
    extrinsics = [AffineMap{Matrix{Float64}, Vector{Float64}}(collect(I(3)), [0, 0, 1])]
    scale = LinearMap{Matrix{Float64}}(scale_factor .* collect(I(3)))
    k = zeros(3)
    files = ["filename"]
    Calibration(CalibrationIO(intrinsic, extrinsics, scale, k, files))
end

# function filter_files(c; improve_n = nothing, improve_threshold = nothing, 
# """
#     improve
# Identify all the images that had relatively high reprojection errors, and rerun the calibration without them. Include a maximum of `n` images with the lowest reprojection error, or all the images with an error lower than `threshold`.
# """
# function improve(cf, n, threshold, aspect, with_distortion)
#     n_files = length(cf.files)
#     n_files ≤ n && return cf
#     reprojection = sqrt.(_reprojection.(Ref(cf), 1:n_files) ./ prod(cf.n_corners))
#     cutoff = max(threshold, sort(reprojection)[n])
#     files = [file for (file, ϵ) in zip(cf.files, reprojection) if ϵ ≤ cutoff]
#     CalibrationFit(files, cf.n_corners, cf.checker_size, aspect, with_distortion)
# end
# improve(cf, n::Int, threshold::Nothing, aspect, with_distortion) = improve(cf, n, 2, aspect, with_distortion)
# improve(cf, n::Nothing, threshold::Int, aspect, with_distortion) = improve(cf, 15, threshold, aspect, with_distortion)
# improve(cf, n::Nothing, threshold::Nothing, aspect, with_distortion) = cf
