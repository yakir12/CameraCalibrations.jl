function obj2img(Rs, ts, frow, fcol, crow, ccol, checker_size)
    intrinsic = AffineMap(SDiagonal(frow, fcol), SVector(crow, ccol))
    extrinsics = AffineMap.(Base.splat(RotationVec).(Rs), SVector{3, Float64}.(ts))
    scale = LinearMap(SDiagonal{3}(I/checker_size))
    return intrinsic, extrinsics, scale
end

"""
    Calibration(files, n_corners, checker_size, extrinsic_index)
Build a calibration object. `files` are the image files of the checkerboard. `n_corners` is a tuple of the number of corners in each of the sides of the checkerboard. `checker_size` is the physical size of the checker (e.g. in cm). `with_distortion` controls if radial lens distortion is included in the model or not.
"""
function fit(files, n_corners, checker_size; aspect = 1, with_distortion = true, inverse_samples = 100, plot_folder::Union{Nothing, String} = nothing)

    files, objpoints, imgpointss, sz, k, Rs, ts, frow, fcol, crow, ccol = detect_fit(unique(files), n_corners, with_distortion, aspect)
    objpoints .*= checker_size
    intrinsic, extrinsics, scale = obj2img(Rs, ts, frow, fcol, crow, ccol, checker_size)
    c = Calibration(intrinsic, extrinsics, scale, k, files)

    # cf = improve(c, improve_n, improve_threshold, aspect, with_distortion)

    ϵ = calculate_errors(c, imgpointss, objpoints, checker_size, sz, files, n_corners, inverse_samples)
    plot(plot_folder, c, imgpointss, n_corners, checker_size, sz)
    return (c, ϵ)
end

function _reprojection(c, i, objpoints, imgpoints)
    reprojected = c.(objpoints, i)
    sum(LinearAlgebra.norm_sqr, reprojected .- imgpoints)
end

"""
    calculate_errors(c)
Calculate reprojection, projection, distance, and inverse errors for the calibration `c`. `distance` measures the mean error of the distance between all adjacent checkerboard corners from the expected `checker_size`. `inverse` measures the mean error of applying the calibration's transformation and its inverse `inverse_samples` times.
"""
function calculate_errors(c, imgpointss, objpoints, checker_size, sz, files, n_corners, inverse_samples=100)
    # @show c.extrinsics
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
            rc = rand(RowCol{Float64}) .* (sz .- 1) .+ 1
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
