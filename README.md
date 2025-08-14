# CameraCalibrations

[![Build Status](https://github.com/yakir12/CameraCalibrations.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/yakir12/CameraCalibrations.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/yakir12/CameraCalibrations.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/yakir12/CameraCalibrations.jl)
[![Aqua](https://raw.githubusercontent.com/JuliaTesting/Aqua.jl/master/badge.svg)](https://github.com/JuliaTesting/Aqua.jl)

This is a package for camera calibration in Julia.

## How to install
```julia
] add CameraCalibrations
```

## How to use
### Using files

First we build the calibration object based on the `files`: the image files of the checkerboard, `n_corners`: the number of inner corners in each of the sides of the checkerboard, and `checker_size`: the physical size of the checker (e.g. in cm).

```julia
using CameraCalibrations
c, ϵ = fit(files, n_corners, checker_size)
```

Then we can use that object to calibrate pixel coordinates and/or images, given the extrinsic image we want to refer to:
```julia
i1 = RowCol(1.2, 3.4) # a cartesian index in image-pixel coordinates
xyz = c(i1, 1) # convert to real-world coordinates of the first image
i2 = c(xyz, 1) # convert back to pixel coordinates
i2 ≈ i1 # true
```

The error term, `ϵ`, includes the reprojection, projection, distance, and inverse errors for the calibration. `distance` measures the mean error of the distance between all adjacent checkerboard corners from the expected `checker_size`. `inverse` measures the mean error of applying the calibration's transformation and its inverse 100 times.

### Using im-memory images
The syntax for in-memory images is very similar to that of files. The difference is that here we must specify "`tags`" (a name in the form of a string) for each of the images. Because corner detection can fail in one or more of the images (due to occlusion, low image quality, etc), it can be important for the user to know *which* of the images failed. Thsi becomes important if say one of the images that failed was the one you had hoped to use for the extrinsic parameters. The images must also be gray-scale.

```julia
fit(tags::Vector{T}, imgs::Vector{Matrix{S}}, n_corners, checker_size) where {T <: AbstractString, S <: Gray}
```

### Additional arguments
- `aspect`: Specifies the aspect ratio of the images. Defaults to `1`.
- `with_distortion`: Include lens distortion in the model. This can be useful to exclude if the resulting camera model results in "donut artifacts" (where the projected coordinates wrap back on themselves at the periphery of the image). Defaults to `true`.
- `plot_folder`: Save the rectified calibration images with a red cross on each detected checkerboard corner and a blue one for the reprojected one. This is useful for assessing the quality of the calibration: the checkerboards should look square and the centers of the red and blue crosses should overlap.

## Features
- [x] thread safe
- [x] saving and loading (JSON) calibration files
- [x] corner detection is done with opencv
- [x] model fitting is done with opencv
- [x] opencv is python-free, via `OpenCV.jl`
- [x] plot calibrated images
- [x] allows for calibration images that were saved with an aspect ration ≠ 1
- [x] in-memory images

## Citing
See [`CITATION.bib`](CITATION.bib) for the relevant reference(s).
