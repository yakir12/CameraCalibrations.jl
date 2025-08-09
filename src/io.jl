StructTypes.StructType(::Type{<:Diagonal}) = StructTypes.CustomStruct()
StructTypes.lower(x::Diagonal) = x.diag
StructTypes.lowertype(::Type{Diagonal{T, V}}) where {T, V} = V

struct CalibrationIO
    intrinsic::AM
    extrinsics::Vector{AMext}
    scale::LM3
    k::Float64
    files::Vector{String}
end

CalibrationIO(c::Calibration) = CalibrationIO(c.intrinsic, c.extrinsics, c.scale, c.k, c.files)

function Calibration(cio::CalibrationIO)
    distort(rc) = lens_distortion(rc, cio.k)
    real2image = .∘(Ref(cio.intrinsic), distort, Ref(PerspectiveMap()), cio.extrinsics, Ref(cio.scale))
    inv_scale, inv_extrinsics, inv_perspective_maps, inv_distort, inv_intrinsic = img2obj(cio.intrinsic, cio.extrinsics, cio.scale, cio.k)
    image2real = .∘(Ref(inv_scale), inv_extrinsics, inv_perspective_maps, inv_distort, Ref(inv_intrinsic))
    return Calibration(cio.intrinsic, cio.extrinsics, cio.scale, cio.k, cio.files, real2image, image2real)
end

Base.rand(rng::AbstractRNG, ::Random.SamplerType{CalibrationIO}) = CalibrationIO(AffineMap(Diagonal(rand(rng, SVector{2, Float64})), rand(rng, SVector{2, Float64})), [AffineMap(rand(rng, RotationVec{Float64}), rand(rng, SVector{3, Float64})) for _ in 1:5], LinearMap(Diagonal(rand(rng, SVector{3, Float64}))), rand(rng), [String(rand(rng, 'a':'z', 5)) for _ in 1:5])
Base.rand(rng::AbstractRNG, ::Random.SamplerType{Calibration}) = Calibration(rand(CalibrationIO))

function load(file)
    cio = JSON3.read(read(file, String), CalibrationIO)
    Calibration(cio)
end

save(file, cio::CalibrationIO) = open(file, "w") do io
    JSON3.write(io, cio)
end

save(file, c::Calibration) = save(file, CalibrationIO(c))
