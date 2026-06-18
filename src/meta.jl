const AM = AffineMap{SDiagonal{2, Float64}, SVector{2, Float64}}
const AMext = AffineMap{RotationVec{Float64}, SVector{3, Float64}}
const LM3 = LinearMap{SDiagonal{3, Float64}}

"""
    RowCol(row, col)
An alias for a static vector of two, row and column, indicating a cartesian coordinate in an image/matrix.
"""
const RowCol = SVector{2, Float32}

"""
    XYZ(x, y, z)
An alias for a static vector of three, x, y, and z, indicating a real-world coordinate. Note that `x` is equivalent to the `column` in `RowCol` and the `y` is equivalent to the `row`.
"""
const XYZ = SVector{3, <: Real}

struct Calibration
    intrinsic::AM
    extrinsics::Vector{AMext}
    scale::LM3
    k::Vector{Float64}
    files::Vector{String}
    real2image::Vector{ComposedFunction}
    image2real::Vector{ComposedFunction}
    center::Union{Nothing, RowCol}
    north::Union{Nothing, RowCol}
    extrinsic_index::Union{Nothing, Int}
end

# Minimum allowed distance, in pixels, between `center` and `north`. Below this the direction `center → north` is too
# short to estimate the alignment angle reliably (a small pixel error implies a large angular error).
const MIN_NORTH_DISTANCE = 10

"""
    _validate_reference(extrinsic_index, center, north, n_extrinsics)
Check that the optional reference-frame metadata is a valid combination. The rule is north ⇒ center ⇒
`extrinsic_index`: a `north` direction requires both a `center` and an `extrinsic_index`, and a `center` requires an
`extrinsic_index`. An `extrinsic_index` alone (or no metadata at all) is fine. When given, `extrinsic_index` must be
a valid index into the `n_extrinsics` extrinsic images, and `north` must be at least `MIN_NORTH_DISTANCE` pixels from
`center`.
"""
function _validate_reference(extrinsic_index, center, north, n_extrinsics)
    if !isnothing(north) && (isnothing(center) || isnothing(extrinsic_index))
        error("`north` requires both `center` and `extrinsic_index`")
    elseif !isnothing(center) && isnothing(extrinsic_index)
        error("`center` requires `extrinsic_index`")
    end
    if !isnothing(extrinsic_index) && !(1 ≤ extrinsic_index ≤ n_extrinsics)
        error("`extrinsic_index` ($extrinsic_index) must be between 1 and the number of extrinsic images ($n_extrinsics)")
    end
    if !isnothing(north) && norm(north - center) < MIN_NORTH_DISTANCE
        error("`north` is only $(round(norm(north - center), digits = 1)) pixels from `center`; they must be at least $MIN_NORTH_DISTANCE pixels apart for a reliable direction")
    end
    return nothing
end

function Calibration(intrinsic, extrinsics, scale, k, files; center = nothing, north = nothing, extrinsic_index = nothing)
    center = isnothing(center) ? nothing : RowCol(center)
    north = isnothing(north) ? nothing : RowCol(north)
    _validate_reference(extrinsic_index, center, north, length(extrinsics))
    distort(rc) = lens_distortion(rc, k)
    real2image = .∘(Ref(intrinsic), distort, Ref(PerspectiveMap()), extrinsics, Ref(scale))
    inv_scale, inv_extrinsics, inv_perspective_maps, inv_distort, inv_intrinsic = img2obj(intrinsic, extrinsics, scale, k)
    image2real = .∘(Ref(inv_scale), inv_extrinsics, inv_perspective_maps, inv_distort, Ref(inv_intrinsic))
    Calibration(intrinsic, extrinsics, scale, k, files, real2image, image2real, center, north, extrinsic_index)
end

"""
    lens_distortion
Lens distortion for up to 3 radial coefficients.
"""
function lens_distortion(v, k)
    kk = zeros(7)
    kk[1] = 1
    for (i, ki) in zip(3:2:7, k)
        kk[i] = ki
    end
    p = Polynomial(kk)
    return v*p(norm(v))
end


"""
    inv_lens_distortion
Inverse lens distortion for up to 3 radial coefficients.
"""
function inv_lens_distortion(v2, k)
    f(v) = norm(lens_distortion(v, k) - v2)
    x0 = [v2...]
    o = optimize(f, x0)
    SVector{2, Float64}(o.minimizer)
end

# this is the inverse prespective map
depth(rc1, t, l) = -t/(l⋅rc1)
function get_inv_prespective_map(inv_extrinsic)
    function (rc)
        rc1 = push(rc, 1)
        t = inv_extrinsic.translation[3]
        l = inv_extrinsic.linear[end, :]
        d = depth(rc1, t, l)
        return d .* rc1
    end
end

function img2obj(intrinsic, extrinsics, scale, k)
    inv_extrinsics = inv.(extrinsics)
    inv_perspective_maps = get_inv_prespective_map.(inv_extrinsics)
    inv_distort(rc) = inv_lens_distortion(rc, k)
    return inv(scale), inv_extrinsics, inv_perspective_maps, inv_distort, inv(intrinsic)
end

# Bounds-check a directly-supplied extrinsic index, with a clear message instead of a bare BoundsError.
function _check_extrinsic_index(c::Calibration, extrinsic_index)
    1 ≤ extrinsic_index ≤ length(c.extrinsics) || throw(ArgumentError("`extrinsic_index` ($extrinsic_index) must be between 1 and the number of extrinsic images ($(length(c.extrinsics)))"))
    return nothing
end

"""
    c(i::RowCol, extrinsic)
Convert the row column StaticArray `i` to its real-world equivalent, `XYZ`, for the extrinsic parameters from the `extrinsic` image (given as an index or file name).
"""
function (c::Calibration)(i::RowCol, extrinsic_index::Int)
    _check_extrinsic_index(c, extrinsic_index)
    c.image2real[extrinsic_index](i)
end

"""
    c(xyz::XYZ, extrinsic)
Convert the x, y, z StaticArray `xyz` to its pixel-coordinate equivalent, `RowCol` for the extrinsic parameters from the `extrinsic` image (given as an index or file name).
"""
function (c::Calibration)(xyz::XYZ, extrinsic_index::Int)
    _check_extrinsic_index(c, extrinsic_index)
    c.real2image[extrinsic_index](xyz)
end

"""
    rectification(c, extrinsic_index)
Return a function that accepts an instance of `::RowCol` and converts it to its real-world equivalent for the extrinsic parameters from the `extrinsic_index` image, without its third dimension (which would be ≈ 0): an `xy` coordinate.
"""
function rectification(c::Calibration, extrinsic_index)
    _check_extrinsic_index(c, extrinsic_index)
    pop ∘ c.image2real[extrinsic_index]
end

"""
    reference_affine(c)
The planar (2D) part of `rectification(c)`: a `Translation` placing `center` at the origin `(0, 0)` and, when a
`north` is set, a rotation aligning the `center → north` direction to `+y` (up). It is the single source of truth for
the reference-frame geometry, so that `rectification(c) == reference_affine(c) ∘ pop ∘ image2real[extrinsic_index]`
and the rectified-image plot share exactly the same transform. Assumes `center` (hence `extrinsic_index`) is set.
Returned as an invertible `CoordinateTransformations` map.
"""
function reference_affine(c::Calibration)
    f = pop ∘ c.image2real[c.extrinsic_index]
    cc = f(c.center)
    A = Translation(-cc)
    isnothing(c.north) && return A
    p = f(c.north) - cc
    return LinearMap(Angle2d(π/2 - atan(p[2], p[1]))) ∘ A
end

"""
    rectification(c)
Return a function that accepts an instance of `::RowCol` and converts it to its real-world `xy` coordinate in the
reference frame defined by the calibration's `center`, `north`, and `extrinsic_index` metadata. Requires at least an
`extrinsic_index`. When a `center` is set it becomes the origin `(0, 0)`; when a `north` is also set the `center →
north` direction is rotated to point towards `+y` (up). With only an `extrinsic_index` this is the plain deskew,
equivalent to `rectification(c, c.extrinsic_index)`.
"""
function rectification(c::Calibration)
    isnothing(c.extrinsic_index) && error("rectification(c) requires an `extrinsic_index`; use rectification(c, i) for a plain deskew")
    base = pop ∘ c.image2real[c.extrinsic_index]
    isnothing(c.center) && return base
    return reference_affine(c) ∘ base
end

