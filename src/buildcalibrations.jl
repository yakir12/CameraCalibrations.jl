function obj2img(Rs, ts, frow, fcol, crow, ccol, checker_size)
    intrinsic = AffineMap(SDiagonal(frow, fcol), SVector(crow, ccol))
    extrinsics = AffineMap.(Base.splat(RotationVec).(Rs), SVector{3, Float64}.(ts))
    scale = LinearMap(SDiagonal{3}(I/checker_size))
    return intrinsic, extrinsics, scale
end

# Identify the chosen extrinsic image by its tag (file name), before detection drops any image. Accepts either an
# integer position into the original `tags` or one of the tags directly. Returns the tag, or `nothing`.
_extrinsic_tag(::Nothing, tags) = nothing
function _extrinsic_tag(extrinsic_index::Integer, tags)
    1 ≤ extrinsic_index ≤ length(tags) || error("`extrinsic_index` ($extrinsic_index) must be between 1 and the number of images ($(length(tags)))")
    return tags[extrinsic_index]
end
function _extrinsic_tag(extrinsic_index::AbstractString, tags)
    extrinsic_index in tags && return extrinsic_index
    # fall back to matching by file name (basename), e.g. when `tags` are full paths
    matches = findall(t -> basename(t) == extrinsic_index, tags)
    isempty(matches) && error("`extrinsic_index` (\"$extrinsic_index\") is not one of the image tags/file names")
    length(matches) > 1 && error("`extrinsic_index` (\"$extrinsic_index\") matches more than one image; pass the full tag to disambiguate")
    return tags[only(matches)]
end

# After detection (which drops images whose corners weren't found), map the chosen extrinsic tag to its position in
# the kept `files`. Errors if corners were not detected in that image, since it cannot serve as the reference.
_resolve_extrinsic_index(::Nothing, files) = nothing
function _resolve_extrinsic_index(extrinsic_tag, files)
    i = findfirst(==(extrinsic_tag), files)
    isnothing(i) && error("corners were not detected in the chosen extrinsic image (\"$extrinsic_tag\"), so it cannot be used as the extrinsic reference")
    return i
end

"""
    Calibration(tags::Vector{AbstractString}, imgs::Vector{Matrix{Gray}}, n_corners, checker_size; aspect = 1, radial_parameters = 1, plot_folder::Union{Nothing, String} = nothing)

Build a calibration object. `tags` are references (just names) to each of the images in `imgs`. Typically, one of these will be the chosen image for the extrinsics parameters. `imgs` are the images of the checkerboard. `n_corners` is a tuple of the number of corners in each of the sides of the checkerboard. `checker_size` is the physical size of the checker (e.g. in cm). `aspect` is the aspect ratio of the images in `imgs` (most commonly 1). `radial_parameters` controls how many radial lens distortion parameters are included in the model. `extrinsic_index` selects the reference image either by its position in `tags` (an `Int`) or by its tag/file name (a `String`); it is resilient to images dropped during corner detection. Finally, `plot_folder` is the path to a directory where diagnostic images can be saved to.
"""
function fit(tags::Vector{T}, imgs::Vector{Matrix{S}}, n_corners, checker_size; aspect = 1, radial_parameters::Int = 1, plot_folder::Union{Nothing, String} = nothing, center = nothing, north = nothing, extrinsic_index = nothing) where {T <: AbstractString, S <: Gray}
    length(tags) == length(imgs) || throw(ArgumentError("`tags` and `imgs` should have the same length"))
    extrinsic_tag = _extrinsic_tag(extrinsic_index, tags)
    files, objpoints, imgpointss, sz, k, Rs, ts, frow, fcol, crow, ccol, extrinsic_index, keep = detect_fit(tags, imgs, n_corners, radial_parameters, aspect, extrinsic_tag)
    kept_imgs = imgs[keep]   # the images whose corners were detected, aligned with `files`/`imgpointss`
    objpoints = objpoints .* checker_size
    intrinsic, extrinsics, scale = obj2img(Rs, ts, frow, fcol, crow, ccol, checker_size)
    c = Calibration(intrinsic, extrinsics, scale, k, files; center, north, extrinsic_index)

    # cf = improve(c, improve_n, improve_threshold, aspect, with_distortion)

    ϵ = calculate_errors(c, imgpointss, objpoints, checker_size, sz, files, n_corners)
    plot(plot_folder, c, kept_imgs, imgpointss, n_corners, checker_size, sz)
    if plot_folder !== nothing && c.center !== nothing && c.extrinsic_index !== nothing
        plot_rectification(plot_folder, c, kept_imgs, imgpointss, n_corners, checker_size)
    end
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

function only_scale(scale_factor; center = nothing, north = nothing, extrinsic_index = nothing)
    intrinsic = AffineMap{Matrix{Float64}, Vector{Float64}}(collect(I(2)), zeros(2))
    extrinsics = [AffineMap{Matrix{Float64}, Vector{Float64}}(collect(I(3)), [0, 0, 1])]
    scale = LinearMap{Matrix{Float64}}(scale_factor .* collect(I(3)))
    k = zeros(3)
    files = ["filename"]
    Calibration(CalibrationIO(intrinsic, extrinsics, scale, k, files; center, north, extrinsic_index))
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
