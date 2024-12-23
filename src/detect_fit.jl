# the following convertion functions are necessary due to some overly eager type checking in the python functions.
convert_from_py_corners(py_corners, n_corners) = reshape(RowCol.(eachslice(PyArray(py_corners); dims = 1)), n_corners)

convert_to_py_objpoints(objpoints, n) = PyList(fill(np.array(Float32.(reshape(reduce((x1, x2) -> hcat(x1, Vector(x2)), objpoints)', (1, length(objpoints), 3)))), n))

convert_to_py_imgpointss(imgpointss) = PyList([np.array(reshape(reduce((x1, x2) -> hcat(x1, Vector(x2)), imgpoints)', 1, length(imgpoints), 2)) for imgpoints in imgpointss])

"""
    get_object_points
Produce the real-world locations of the corners of the checkerboard.
"""
function get_object_points(n_corners)
    objpoints = Matrix{XYZ}(undef, n_corners)
    for i in CartesianIndices(n_corners)
        x, y = Tuple(i) .- 1
        objpoints[i] = XYZ(x, y, 0)
    end
    return objpoints
end

"""
    _detect_corners
Wraps OpenCV function to auto-detect corners in an image.
"""
function _detect_corners(file, n_corners)
    criteria = (cv2.TERM_CRITERIA_EPS + cv2.TERM_CRITERIA_MAX_ITER, 30, 0.001)

    img = cv2.imread(file)
    gray = cv2.cvtColor(img, cv2.COLOR_BGR2GRAY)
    ret, py_corners = cv2.findChessboardCorners(gray, n_corners, cv2.CALIB_CB_ADAPTIVE_THRESH + cv2.CALIB_CB_NORMALIZE_IMAGE + cv2.CALIB_CB_FAST_CHECK)

    !Bool(ret) && return missing

    ref_corners = cv2.cornerSubPix(gray, py_corners, (11,11),(-1,-1), criteria)
    flp_corners = np.flip(ref_corners, axis = 2)
    corners = convert_from_py_corners(flp_corners, n_corners)
    return (file, reverse(corners, dims = 1))
end

function detect_corners(_files, n_corners)
    fi = skipmissing(_detect_corners.(_files, Ref(n_corners)))
    return (; files = first.(fi), imgpointss = last.(fi))
end

"""
    fit_model
Wraps OpenCV function to fit a camera model to given object and image points.
"""
function fit_model(sz, objpoints, imgpointss, n_corners,  with_distortion, aspect)
    flags = cv2.CALIB_ZERO_TANGENT_DIST + cv2.CALIB_FIX_K3 + cv2.CALIB_FIX_K2 + (with_distortion ? 0 : cv2.CALIB_FIX_K1) + cv2.CALIB_FIX_ASPECT_RATIO

    cammatrix = convert(Matrix{Float32}, I(3))
    cammatrix[1] = aspect

    cammatrix = np.ascontiguousarray(np.array(cammatrix)) # https://stackoverflow.com/a/50128836/2261957

    _, py_mtx, py_dist, py_rvecs, py_tvecs = cv2.calibrateCamera(convert_to_py_objpoints(objpoints, length(imgpointss)), convert_to_py_imgpointss(imgpointss), np.flip(sz), cammatrix, nothing; flags)

    k, _ = PyArray(py_dist)
    @assert with_distortion || k == 0 "distortion was $with_distortion but k isn't zero:" k

    Rs = [vec(pyconvert(Matrix{Float64}, x)) for x in py_rvecs]

    ts = [vec(pyconvert(Matrix{Float64}, x)) for x in py_tvecs]

    mtx = Matrix(PyArray(py_mtx))
    frow = mtx[1,1]
    fcol = mtx[2,2]
    crow = mtx[1,3]
    ccol = mtx[2,3]

    return (; k, Rs, ts, frow, fcol, crow, ccol)
end

function detect_fit(_files, n_corners, with_distortion, aspect)
    files, imgpointss = detect_corners(_files, n_corners)

    objpoints = get_object_points(n_corners)
    sz = size(FileIO.load(files[1]))

    k, Rs, ts, frow, fcol, crow, ccol = fit_model(sz, objpoints, imgpointss, n_corners, with_distortion, aspect)

    return (; files, objpoints, imgpointss, sz, k, Rs, ts, frow, fcol, crow, ccol)
end



