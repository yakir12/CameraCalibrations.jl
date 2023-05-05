"""
    _detect_corners
Wraps OpenCV function to auto-detect corners in an image.
"""
function _detect_corners(file, n_corners)
    img = load(file)
    gray = np.array(rawview(channelview(img)))
    ret, py_corners = cv2.findChessboardCorners(gray, n_corners, cv2.CALIB_CB_ADAPTIVE_THRESH + cv2.CALIB_CB_NORMALIZE_IMAGE + cv2.CALIB_CB_FAST_CHECK)
    !Bool(ret) && return missing
    ref_corners = cv2.cornerSubPix(gray, py_corners, (11,11),(-1,-1), criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 30, 0.001))
    corners = reshape(reverse.(SV.(eachslice(PyArray(ref_corners); dims=1))), n_corners) # reversed cause these are not cartesian indices
    return (file, corners)
end

function detect_corners(_files, n_corners)
    fi = skipmissing(_detect_corners.(_files, Ref(n_corners)))
    @assert !isempty(fi) "No corners were detected in any of the image files"
    return (first.(fi), last.(fi))
end

"""
    fit_model
Wraps OpenCV function to fit a camera model to given object and image points.
"""

function Calibration(_files, n_corners, checker_size; lens_distortion = LensDistortion[])
    # if !issorted(n_corners)
    #     n_corners = reverse(n_corners)
    # end
    files, images_points = detect_corners(_files, n_corners)

    object_points = [reverse(checker_size .* SV(Tuple(ij))) .- 1 for ij in CartesianIndices(n_corners)] # reversed to keep the warped images not mirrored

    Calibration(files, n_corners, checker_size, images_points, object_points, size(load(files[1])), lens_distortion)
end

function Camera(calib::Calibration)
    flags = mapreduce(k -> lens_distortions[k], +, calib.distortion, init=0)
    # flags = cv2.CALIB_ZERO_TANGENT_DIST + cv2.CALIB_FIX_K3 + cv2.CALIB_FIX_K2 + (calib.with_distortion ? 0 : cv2.CALIB_FIX_K1)

    objpoints = fill(np.array(reduce((v, x) -> vcat(v, Float32[x... 0]), calib.object_points ./ calib.checker_size, init=Matrix{Float32}(undef, 0, 3))), length(calib.files))
    imgpoints = PyList(np.array.(reduce.((v, x) -> cat(v, reshape(Float32.(x), 1, 1, 2); dims=1), calib.images_points, init = Array{Float32, 3}(undef, 0, 1, 2))))

    _, py_mtx, py_dist, py_rvecs, py_tvecs = cv2.calibrateCamera(objpoints, imgpoints, calib.image_size, nothing, nothing; flags)

    k₁, k₂, p₁, p₂, k₃ = PyArray(py_dist)
    # @assert calib.with_distortion || k == 0 "distortion was $(calib.with_distortion) but k isn't zero:" k

    Rs = SVector{3, Float64}.(Matrix.(PyArray.(py_rvecs)))

    ts = SVector{3, Float64}.(Matrix.(PyArray.(py_tvecs)))

    mtx = PyMatrix(py_mtx)
    focal_length = (mtx[1,1], mtx[2,2])
    principal_point = SVector{2, Float64}(mtx[1,3], mtx[2,3])

    Camera(1/calib.checker_size, (; k₁, k₂, k₃, p₁, p₂), Rs, ts, focal_length, principal_point)
end

