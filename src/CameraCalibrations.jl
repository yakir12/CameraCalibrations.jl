module CameraCalibrations

using LinearAlgebra, Statistics, Random
import ImageIO, FileIO
using PythonCall, ImageBase, StaticArrays, ImageDraw, ImageTransformations
using Rotations, CoordinateTransformations, Polynomials
using OpenCV_jll

using JSON3, StructTypes

export fit, Calibration, CalibrationIO, rectification, RowCol, XYZ

include("meta.jl")
include("io.jl")
include("detect_fit.jl")
include("buildcalibrations.jl")
include("plot_calibration.jl")

# using CondaPkg
# CondaPkg.add.(["numpy", "opencv"])

const cv2 = PythonCall.pynew()
const np = PythonCall.pynew()

function __init__()
    PythonCall.pycopy!(cv2, pyimport("cv2"))
    PythonCall.pycopy!(np, pyimport("numpy"))
end

end
