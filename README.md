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
First we build the calibration object based on the `files`: the image files of the checkerboard, `n_corners`: the number of inner corners in each of the sides of the checkerboard, and `checker_size`: the physical size of the checker (e.g. in cm).

```julia
using CameraCalibrations
c = fit(files, n_corners, checker_size)
```

Then we can use that object to calibrate pixel coordinates and/or images, given the extrinsic image we want to refer to:
```julia
i1 = RowCol(1.2, 3.4) # a cartesian index in image-pixel coordinates
xyz = c(i1, 1) # convert to real-world coordinates of the first image
i2 = c(xyz, 1) # convert back to pixel coordinates
i2 ≈ i1 # true
```

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
