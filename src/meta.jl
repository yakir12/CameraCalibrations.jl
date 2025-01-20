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
    k::Float64
    files::Vector{String}
    real2image::Vector{ComposedFunction}
    image2real::Vector{ComposedFunction}
end

function Calibration(intrinsic, extrinsics, scale, k, files)
    distort(rc) = lens_distortion(rc, k)
    real2image = .∘(Ref(intrinsic), distort, Ref(PerspectiveMap()), extrinsics, Ref(scale))
    inv_scale, inv_extrinsics, inv_perspective_maps, inv_distort, inv_intrinsic = img2obj(intrinsic, extrinsics, scale, k)
    image2real = .∘(Ref(inv_scale), inv_extrinsics, inv_perspective_maps, inv_distort, Ref(inv_intrinsic))
    Calibration(intrinsic, extrinsics, scale, k, files, real2image, image2real)
end

"""
    lens_distortion
Lens distortion for one radial coefficient.
"""
function lens_distortion(v, k)
    k == 0 && return v
    r² = LinearAlgebra.norm_sqr(v)
    radial = 1 + k*r²
    radial*v
end

"""
    inv_lens_distortion
Analytical inverse lens distortion for one radial coefficient.
"""
function inv_lens_distortion(v2, k)
    k == 0 && return v2
    c = k*LinearAlgebra.norm_sqr(v2)
    rs = roots(Polynomial([c, 0, 1, -1]))
    rrs = filter(x -> abs(imag(x)) < 1e-10  , rs)
    radial = maximum(real, rrs)
    v2 / radial
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

"""
    c(i::RowCol, extrinsic)
Convert the row column StaticArray `i` to its real-world equivalent, `XYZ`, for the extrinsic parameters from the `extrinsic` image (given as an index or file name).
"""
(c::Calibration)(i::RowCol, extrinsic_index::Int) = c.image2real[extrinsic_index](i)

"""
    c(xyz::XYZ, extrinsic)
Convert the x, y, z StaticArray `xyz` to its pixel-coordinate equivalent, `RowCol` for the extrinsic parameters from the `extrinsic` image (given as an index or file name).
"""
(c::Calibration)(xyz::XYZ, extrinsic_index::Int) = c.real2image[extrinsic_index](xyz)

function (c::Calibration)(xyz_or_i)
    extrinsic_index = findfirst(contains(r"extrinsic"), c.files)
    return c(xyz_or_i, extrinsic_index)
end

"""
    rectification(c, extrinsic_index)
Return a function that accepts an instance of `::RowCol` and converts it to its real-world equivalent for the extrinsic parameters from the `extrinsic_index` image, without its third dimension (which would be ≈ 0): an `xy` coordinate.
"""
rectification(c::Calibration, extrinsic_index) = pop ∘ c.image2real[extrinsic_index]
function rectification(c::Calibration) 
    extrinsic_index = findfirst(contains(r"extrinsic"), c.files)
    return pop ∘ c.image2real[extrinsic_index]
end

