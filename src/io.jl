struct CalibrationIO
    intrinsic::AffineMap{Matrix{Float64}, Vector{Float64}}
    extrinsics::Vector{AffineMap{Matrix{Float64}, Vector{Float64}}}
    scale::LinearMap{Matrix{Float64}}
    k::Vector{Float64}
    files::Vector{String}
    center::Union{Nothing, Vector{Float64}}
    north::Union{Nothing, Vector{Float64}}
    extrinsic_index::Union{Nothing, Int}
end

# Convenience constructor so existing positional calls keep working and gain the reference-frame kwargs. Coexists
# with the auto-generated 8-argument constructor used by JSON's typed parsing.
CalibrationIO(intrinsic, extrinsics, scale, k, files; center = nothing, north = nothing, extrinsic_index = nothing) =
    CalibrationIO(intrinsic, extrinsics, scale, k, files, center, north, extrinsic_index)

function convert_intrinsic_or_extrinsic(x::AffineMap{Matrix{T}, Vector{T}}) where {T}
    if length(x.translation) == 2
        AffineMap(SDiagonal(SVector{2, T}(diag(x.linear))), SVector{2, T}(x.translation))
    else
        AffineMap(RotationVec{T}(x.linear), SVector{3, T}(x.translation))
    end
end

convert_intrinsic_or_extrinsic(x::AffineMap{R, SVector{N, T}}) where {R, N, T} = AffineMap(Matrix{T}(x.linear), Vector{T}(x.translation))

convert_scale(scale::LinearMap{Matrix{T}}) where {T} = LinearMap(SDiagonal(SVector{3, T}(diag(scale.linear))))
convert_scale(scale::LinearMap{Diagonal{T, SVector{N, T}}}) where {T, N} = LinearMap(Matrix{T}(scale.linear))

CalibrationIO(c::Calibration) = CalibrationIO(convert_intrinsic_or_extrinsic(c.intrinsic), convert_intrinsic_or_extrinsic.(c.extrinsics), convert_scale(c.scale), c.k, c.files, isnothing(c.center) ? nothing : Vector{Float64}(c.center), isnothing(c.north) ? nothing : Vector{Float64}(c.north), c.extrinsic_index)

function Calibration(cio::CalibrationIO)
    distort(rc) = lens_distortion(rc, cio.k)
    intrinsic = convert_intrinsic_or_extrinsic(cio.intrinsic)
    scale = convert_scale(cio.scale)
    extrinsics = convert_intrinsic_or_extrinsic.(cio.extrinsics)
    real2image = .∘(Ref(intrinsic), distort, Ref(PerspectiveMap()), extrinsics, Ref(scale))
    inv_scale, inv_extrinsics, inv_perspective_maps, inv_distort, inv_intrinsic = img2obj(intrinsic, extrinsics, scale, cio.k)
    image2real = .∘(Ref(inv_scale), inv_extrinsics, inv_perspective_maps, inv_distort, Ref(inv_intrinsic))
    center = isnothing(cio.center) ? nothing : RowCol(cio.center)
    north = isnothing(cio.north) ? nothing : RowCol(cio.north)
    _validate_reference(cio.extrinsic_index, center, north, length(extrinsics))
    return Calibration(intrinsic, extrinsics, scale, cio.k, cio.files, real2image, image2real, center, north, cio.extrinsic_index)
end

Base.rand(rng::AbstractRNG, ::Random.SamplerType{CalibrationIO}) = CalibrationIO(AffineMap(rand(rng, 2, 2), rand(rng, 2)), [AffineMap(rand(rng, 3, 3), rand(rng, 3)) for _ in 1:5], LinearMap(rand(rng, 3, 3)), rand(rng, 3), [String(rand(rng, 'a':'z', 5)) for _ in 1:5]; center = rand(rng, 2), north = rand(rng, 2), extrinsic_index = rand(rng, 1:5))
Base.rand(rng::AbstractRNG, ::Random.SamplerType{Calibration}) = Calibration(rand(rng, CalibrationIO))

_ismatlab(file) = read(file, 6) == codeunits("MATLAB")

"""
    load(file)
Load a `Calibration` from a JSON file written by [`save`](@ref). All reference-frame metadata (`center`, `north`,
`extrinsic_index`) is read from the file. To load a MATLAB Camera Calibrator `.mat` file, use [`loadMAT`](@ref)
instead.
"""
function load(file)
    _ismatlab(file) && error("\"$file\" looks like a MATLAB file; use `loadMAT` to load MATLAB calibrations")
    Calibration(JSON.parsefile(file, CalibrationIO))
end

save(file, cio::CalibrationIO) = JSON.json(file, cio)

save(file, c::Calibration) = JSON.json(file, CalibrationIO(c))


"""
    loadMAT(file; extrinsic_index = 1, center = nothing, north = nothing)
Load a `Calibration` from a MATLAB Camera Calibrator `.mat` file (read-only). Since `.mat` files do not carry this
package's reference-frame metadata, supply it here: `extrinsic_index` selects the reference image and `center`/`north`
define the reference frame. To load this package's own JSON files, use [`load`](@ref) instead.
"""
function loadMAT(file; extrinsic_index::Int = 1, center = nothing, north = nothing)
    _ismatlab(file) || error("\"$file\" does not look like a MATLAB file; use `load` for this package's JSON calibrations")
    dict = matread(file)
    k = only(keys(dict))
    dict = dict[k]
    for k in ("TranslationVectors", "RotationVectors", "RadialDistortion", "K")
        @assert haskey(dict, k) "MATLAB Camera calibration file missing $k"
    end

    fcol = dict["K"][1,1]
    frow = dict["K"][2,2]
    ccol = dict["K"][1,3]
    crow = dict["K"][2,3]
    intrinsic = AffineMap(Matrix(Diagonal([frow, fcol])), [crow, ccol])

    # both of these have their x and y the other way around, due ot some matlab convension
    Rs = [Matrix{Float64}(RotationVec(-R[[2,1,3]]...)) for R in eachrow(dict["RotationVectors"])] # negative due to some matlab angle convension...
    ts = [Vector{Float64}(t[[2,1,3]]) for t in eachrow(dict["TranslationVectors"])]
    extrinsics = AffineMap.(Rs, ts)

    scale = LinearMap(Matrix(Diagonal{Float64}(I(3))))

    k = vec(dict["RadialDistortion"])

    files =  string.("image ", 1:length(Rs))
    files[extrinsic_index] = "extrinsic"

    cio = CalibrationIO(intrinsic, extrinsics, scale, k, files; center, north, extrinsic_index)
    Calibration(cio)
end
