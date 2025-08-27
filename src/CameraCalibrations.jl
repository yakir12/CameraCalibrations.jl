module CameraCalibrations

using CoordinateTransformations: CoordinateTransformations, AffineMap, LinearMap, PerspectiveMap, ∘
using ImageBase: ImageBase, @colorant_str, Gray, RGB, channelview, norm, rawview, ⋅
using ImageDraw: ImageDraw, Cross, Point, draw!
using ImageTransformations: ImageTransformations, warp
using JSON3: JSON3
using LinearAlgebra: LinearAlgebra, /, Diagonal, I, convert
using OhMyThreads: OhMyThreads, tcollect
using OpenCV: OpenCV
using Polynomials: Polynomials, Polynomial, roots
using Random: Random, AbstractRNG
using Rotations: Rotations, RotationVec
using StaticArrays: StaticArrays, SDiagonal, SVector, pop, push
using Statistics: Statistics, mean
using StructTypes: StructTypes
using FileIO: FileIO
using MAT: matread
using Optim:optimize

const CRITERIA = OpenCV.TermCriteria(OpenCV.TERM_CRITERIA_EPS + OpenCV.TERM_CRITERIA_MAX_ITER, 30, 0.001)

export fit, Calibration, CalibrationIO, rectification, RowCol, XYZ

include("meta.jl")
include("io.jl")
include("detect_fit.jl")
include("buildcalibrations.jl")
include("plot_calibration.jl")

end
