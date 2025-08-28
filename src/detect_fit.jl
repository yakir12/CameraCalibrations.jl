"""
    _detect_corners
Wraps OpenCV function to auto-detect corners in an image.
"""
function _detect_corners(img, n_corners)
    gry = OpenCV.Mat(reshape(rawview(channelview(img)), 1, size(img)...))
    corners = Matrix{RowCol}(undef, n_corners)
    ret, _ = OpenCV.findChessboardCorners(gry, OpenCV.Size{Int32}(n_corners...), OpenCV.Mat(reshape(reinterpret(Float32, corners), 2, 1, prod(n_corners))), OpenCV.CALIB_CB_ADAPTIVE_THRESH + OpenCV.CALIB_CB_NORMALIZE_IMAGE + OpenCV.CALIB_CB_FAST_CHECK)
    return ret ? corners : missing
    # ref_corners = OpenCV.cornerSubPix(gry, cv_corners, OpenCV.Size{Int32}(11,11), OpenCV.Size{Int32}(-1,-1), CRITERIA)
    # corners = reshape(RowCol.(eachslice(ref_corners, dims = 3)), n_corners)
end

"""
    fit_model
Wraps OpenCV function to fit a camera model to given object and image points.
"""
function fit_model(sz, objpoints, imgpointss, n_corners, radial_parameters, aspect)
    cammat = convert(Matrix{Float64}, I(3))
    cammat[1] = 1
    cammat[2, 2] = aspect
    dist = Vector{Float64}(undef, 5)
    nfiles = length(imgpointss)
    r = [Vector{Float64}(undef, 3) for _ in 1:nfiles]
    t = [Vector{Float64}(undef, 3) for _ in 1:nfiles]
    CALIB_FIX_K = sum([OpenCV.CALIB_FIX_K1, OpenCV.CALIB_FIX_K2, OpenCV.CALIB_FIX_K3][setdiff(1:3, 1:radial_parameters)])
    # @show Int(CALIB_FIX_K)
    flags = OpenCV.CALIB_ZERO_TANGENT_DIST + CALIB_FIX_K + OpenCV.CALIB_FIX_ASPECT_RATIO
    criteria = OpenCV.TermCriteria(OpenCV.TERM_CRITERIA_EPS + OpenCV.TERM_CRITERIA_MAX_ITER, 30, 0.001)

    OpenCV.calibrateCamera(OpenCV.InputArray[Float32.(reshape(stack(objpoints), 3, 1, :)) for _ in 1:nfiles], 
                                                         OpenCV.InputArray[Float32.(reshape(stack(imgpoints), 2, 1, :)) for imgpoints in imgpointss], 
                                                         OpenCV.Size{Int32}(sz...),  
                                                         OpenCV.Mat(reshape(cammat, 1, 3, 3)), 
                                                         OpenCV.Mat(reshape(dist, 1, 1, 5)), 
                                                         OpenCV.InputArray[OpenCV.Mat(reshape(ri, 1, 1, 3)) for ri in r], 
                                                         OpenCV.InputArray[OpenCV.Mat(reshape(ti, 1, 1, 3)) for ti in t], flags, CRITERIA)
    return (k = dist[[1,2,5]], Rs = r, ts = t, frow = cammat[1,1], fcol = cammat[2,2], crow = cammat[3,1], ccol = cammat[3,2])
end

function detect_fit(tags, imgs, n_corners, radial_parameters, aspect)
    fi = tcollect(_detect_corners(img, n_corners) for img in imgs)
    keep = findall(!ismissing, fi)
    tags = tags[keep]
    imgpointss = fi[keep]
    @assert !isempty(tags) "No checkers were detected in any of the images, perhaps try a different `n_corners` (e.g. not $n_corners)"
    objpoints = XYZ.(Tuple.(CartesianIndices((0:(n_corners[1] - 1), 0:(n_corners[2] - 1), 0:0))))
    sz = size(imgs[1])
    k, Rs, ts, frow, fcol, crow, ccol = fit_model(sz, objpoints, imgpointss, n_corners, radial_parameters, aspect)
    return (; tags, objpoints, imgpointss, sz, k, Rs, ts, frow, fcol, crow, ccol)
end
