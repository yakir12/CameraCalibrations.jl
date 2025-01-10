module CameraCalibrations

using LinearAlgebra, Statistics, Random
import ImageIO, FileIO
using ImageBase, StaticArrays, ImageDraw, ImageTransformations
using Rotations, CoordinateTransformations, Polynomials
using OpenCV
using JSON3, StructTypes
const CRITERIA = OpenCV.TermCriteria(OpenCV.TERM_CRITERIA_EPS + OpenCV.TERM_CRITERIA_MAX_ITER, 30, 0.001)

export fit, Calibration, CalibrationIO, rectification, RowCol, XYZ

include("meta.jl")
include("io.jl")
include("detect_fit.jl")
include("buildcalibrations.jl")
include("plot_calibration.jl")

end
