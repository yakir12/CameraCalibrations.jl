module CameraCalibrations

using CoordinateTransformations: CoordinateTransformations, AffineMap, LinearMap, PerspectiveMap, ∘
using ImageBase: ImageBase, @colorant_str, Gray, RGB, channelview, norm, rawview, ⋅
using ImageDraw: ImageDraw, Cross, Point, draw!
using ImageTransformations: ImageTransformations, warp
using JSON: JSON
using LinearAlgebra: LinearAlgebra, /, Diagonal, I, convert, diag
using OhMyThreads: OhMyThreads, tcollect
using OpenCV: OpenCV
using Polynomials: Polynomials, Polynomial
using Random: Random, AbstractRNG
using Rotations: Rotations, RotationVec
using StaticArrays: StaticArrays, SDiagonal, SVector, pop, push
using Statistics: Statistics, mean
using FileIO: FileIO
using MAT: matread
using Optim:optimize
import Base: one

const CRITERIA = OpenCV.TermCriteria(OpenCV.TERM_CRITERIA_EPS + OpenCV.TERM_CRITERIA_MAX_ITER, 30, 0.001)

export fit, Calibration, CalibrationIO, rectification, RowCol, XYZ, only_scale

include("meta.jl")
include("io.jl")
include("detect_fit.jl")
include("buildcalibrations.jl")
include("plot_calibration.jl")

end
