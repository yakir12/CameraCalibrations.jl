# StructTypes.StructType(::Type{<:SVector}) = StructTypes.CustomStruct()
# StructTypes.lower(x::SVector) = x.data
# StructTypes.lowertype(::Type{SVector{N, T}}) where {N, T} = NTuple{N, T}

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

#
#
# function load(file)
#
#     file = "tmp.calib"
#
#     open(file, "w") do o
#         JSON.print(o, c)
#     end
#
#     dict = JSON.parsefile(file)
#
#




# using LinearAlgebra, Statistics, Random, StaticArrays, Rotations, CoordinateTransformations, JSON3, StructTypes
#
#
# const AM = AffineMap{SDiagonal{2, Float64}, SVector{2, Float64}}
# const AMext = AffineMap{RotationVec{Float64}, SVector{3, Float64}}
# const LM3 = LinearMap{SDiagonal{3, Float64}}
#
# StructTypes.StructType(::Type{<:Diagonal}) = StructTypes.CustomStruct()
# StructTypes.lower(x::Diagonal) = x.diag
# StructTypes.lowertype(::Type{Diagonal{T, V}}) where {T, V} = V
#
# struct CalibrationIO
#     intrinsic::AM
#     extrinsics::Vector{AMext}
#     scale::LM3
#     k::Float64
#     files::Vector{String}
# end
#
# Random.rand(rng::AbstractRNG, ::Random.SamplerType{CalibrationIO}) = CalibrationIO(AffineMap(Diagonal(rand(rng, SVector{2, Float64})), rand(rng, SVector{2, Float64})), [AffineMap(rand(rng, RotationVec{Float64}), rand(rng, SVector{3, Float64})) for _ in 1:5], LinearMap(Diagonal(rand(rng, SVector{3, Float64}))), rand(rng), [String(rand(rng, 'a':'z', 5)) for _ in 1:5])
#
# x = rand(CalibrationIO)
# v = x.intrinsic
# txt = JSON3.write(v)
#
# txt = "{\"linear\":[0.8369119437634303,0.5233244467702733],\"translation\":[0.4655095816308332,0.6672241801983796]}"
# y = JSON3.read(txt, typeof(v))
#
# using StaticArrays, JSON3, LinearAlgebra
# x = Diagonal(SVector(1, 2))
# txt = JSON3.write(x)
# y = JSON3.read(txt, typeof(x))
#


# module PkgA
# using StructTypes, LinearAlgebra
# StructTypes.StructType(::Type{<:Diagonal}) = StructTypes.CustomStruct()
# StructTypes.lower(x::Diagonal) = x.diag
# StructTypes.lowertype(::Type{Diagonal{T, V}}) where {T, V} = V
# end
#
# using .PkgA, Aqua
# Aqua.test_piracies(PkgA)
#
