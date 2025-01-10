"""
    _detect_corners
Wraps OpenCV function to auto-detect corners in an image.
"""
function _detect_corners(img_file, n_corners)
    gry = OpenCV.imread(img_file, OpenCV.IMREAD_GRAYSCALE)
    corners = Matrix{RowCol}(undef, n_corners)
    ret, _ = OpenCV.findChessboardCorners(gry, OpenCV.Size{Int32}(n_corners...), OpenCV.Mat(reshape(reinterpret(Float32, corners), 2, 1, prod(n_corners))), 0)
    return ret ? (img_file, corners) : missing
    # ref_corners = OpenCV.cornerSubPix(gry, cv_corners, OpenCV.Size{Int32}(11,11), OpenCV.Size{Int32}(-1,-1), CRITERIA)
    # corners = reshape(RowCol.(eachslice(ref_corners, dims = 3)), n_corners)
end

"""
    fit_model
Wraps OpenCV function to fit a camera model to given object and image points.
"""
function fit_model(sz, objpoints, imgpointss, n_corners, with_distortion, aspect)
    cammat = convert(Matrix{Float64}, I(3))
    cammat[1] = aspect
    dist = Vector{Float64}(undef, 5)
    nfiles = length(imgpointss)
    r = [Vector{Float64}(undef, 3) for _ in 1:nfiles]
    t = [Vector{Float64}(undef, 3) for _ in 1:nfiles]
    flags = OpenCV.CALIB_ZERO_TANGENT_DIST + OpenCV.CALIB_FIX_K3 + OpenCV.CALIB_FIX_K2 + (with_distortion ? 0 : OpenCV.CALIB_FIX_K1) + OpenCV.CALIB_FIX_ASPECT_RATIO
    criteria = OpenCV.TermCriteria(OpenCV.TERM_CRITERIA_EPS + OpenCV.TERM_CRITERIA_MAX_ITER, 30, 0.001)

    rms, _cammat, _dist, _r, _t = OpenCV.calibrateCamera(OpenCV.InputArray[Float32.(reshape(stack(objpoints), 3, 1, :)) for _ in 1:nfiles], 
                                                         OpenCV.InputArray[Float32.(reshape(stack(imgpoints), 2, 1, :)) for imgpoints in imgpointss], 
                                                         OpenCV.Size{Int32}(sz...),  
                                                         OpenCV.Mat(reshape(cammat, 1, 3, 3)), 
                                                         OpenCV.Mat(reshape(dist, 1, 1, 5)), 
                                                         OpenCV.InputArray[OpenCV.Mat(reshape(ri, 1, 1, 3)) for ri in r], 
                                                         OpenCV.InputArray[OpenCV.Mat(reshape(ti, 1, 1, 3)) for ti in t], flags, CRITERIA)

    return (; k = dist[1], Rs = r, ts = t, frow = cammat[1,1], fcol = cammat[2,2], crow = cammat[1,3], ccol = cammat[2,3])
end

function detect_fit(_files, n_corners, with_distortion, aspect)
    sz = size(FileIO.load(first(_files)))
    fi = skipmissing(_detect_corners.(_files, Ref(n_corners)))
    @assert !isempty(fi) "No checkers were detected in any of the images, perhaps try a different `n_corners`."
    files = first.(fi)
    imgpointss = last.(fi)
    objpoints = XYZ.(Tuple.(CartesianIndices((0:(n_corners[1] - 1), 0:(n_corners[2] - 1), 0:0))))
    k, Rs, ts, frow, fcol, crow, ccol = fit_model(sz, objpoints, imgpointss, n_corners, with_distortion, aspect)
    return (; files, objpoints, imgpointss, sz, k, Rs, ts, frow, fcol, crow, ccol)
end
