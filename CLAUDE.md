# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

CameraCalibrations.jl is a Julia package for camera calibration using checkerboard images. It wraps OpenCV's calibration functions (via OpenCV.jl) for corner detection and model fitting, then provides bidirectional coordinate transformations between pixel coordinates (RowCol) and real-world coordinates (XYZ).

## Development Commands

### Testing
```bash
# Run all tests (single-threaded)
julia --project -e 'using Pkg; Pkg.test()'

# Run tests with multiple threads (recommended, as package includes thread-safety tests)
julia --project -t 4 -e 'using Pkg; Pkg.test()'

# Run a specific test set
julia --project -t 4 test/runtests.jl
```

### Package Management
```bash
# Install/update dependencies
julia --project -e 'using Pkg; Pkg.instantiate()'

# Build package
julia --project -e 'using Pkg; Pkg.build()'
```

### REPL Development
```julia
# Enter package mode and activate
] activate .
] instantiate

# Load package for interactive testing
using CameraCalibrations
```

## Architecture

### Core Data Flow

1. **Corner Detection** (`detect_fit.jl`): Uses OpenCV to detect checkerboard corners in calibration images
2. **Model Fitting** (`detect_fit.jl`): Fits camera intrinsic/extrinsic parameters using OpenCV's calibrateCamera
3. **Calibration Object** (`meta.jl`, `buildcalibrations.jl`): Constructs bidirectional transformation pipelines
4. **Coordinate Transformations** (`meta.jl`): Applies composed transformations for pixel ↔ real-world conversion

### Key Types and Their Relationships

**RowCol** and **XYZ**: Type aliases for `SVector{2, Float32}` (pixel coordinates) and `SVector{3, Real}` (real-world coordinates).

**Calibration**: The main calibration object containing:
- `intrinsic`: Camera intrinsic parameters (focal length, principal point)
- `extrinsics`: One per calibration image (rotation + translation)
- `scale`: Physical scaling factor from checker_size
- `k`: Radial lens distortion coefficients (up to 3)
- `real2image` and `image2real`: Pre-composed transformation pipelines

**CalibrationIO**: Serialization-friendly version using Matrix/Vector types instead of StaticArrays/RotationVec. Conversion functions in `io.jl` handle bidirectional conversion.

### Transformation Pipeline

**Real-world to pixel** (real2image):
```
XYZ → scale → extrinsic (rotate+translate) → perspective → lens_distortion → intrinsic → RowCol
```

**Pixel to real-world** (image2real):
```
RowCol → inv(intrinsic) → inv_lens_distortion → inv_perspective → inv(extrinsic) → inv(scale) → XYZ
```

Each transformation is a composed function (using `∘`) stored per extrinsic image. The inverse perspective map requires solving for depth using the z-component of the transformed coordinates.

### Thread Safety

The package is thread-safe. Corner detection uses `OhMyThreads.tcollect` for parallel processing. Tests verify thread safety by running calibrations concurrently (`Threads.@threads`).

### Serialization (io.jl)

Supports two formats:
- **JSON**: Uses CalibrationIO for serialization, auto-detected on load
- **MATLAB**: Loads `.mat` files from MATLAB's Camera Calibrator app (read-only)

The CalibrationIO struct uses regular matrices/vectors instead of StaticArrays/RotationVec types for JSON compatibility. Conversion functions handle the mapping between runtime-efficient types (Calibration) and serialization-friendly types (CalibrationIO).

### Lens Distortion

**Forward distortion** (`lens_distortion`): Uses a polynomial with odd powers (1, r³, r⁵, r⁷) to model radial distortion. The `k` vector contains up to 3 radial coefficients.

**Inverse distortion** (`inv_lens_distortion`): Numerically inverts distortion using `Optim.optimize` since there's no closed-form inverse for the polynomial.

### Plotting (plot_calibration.jl)

When `plot_folder` is specified, saves rectified calibration images with:
- Red crosses: detected checkerboard corners
- Blue crosses: reprojected corners from the calibration model

This visualizes calibration quality—crosses should overlap and checkerboards should appear square.

## Important Constraints

### Julia Version Restriction
The package supports Julia 1.0–1.10 only. This is enforced in Project.toml line 47. If updating dependencies or making changes that might affect compatibility, verify against this constraint.

### Aspect Ratio Handling
The `aspect` parameter accounts for images saved with non-square pixels (aspect ≠ 1). This is set during calibration and affects the intrinsic camera matrix (line 21 in detect_fit.jl). If working with aspect ratio code, note that it only affects the y-axis focal length (`cammat[2, 2]`).

### Radial Distortion Parameters
The `radial_parameters` argument (1–3) controls which radial distortion coefficients are fitted vs. fixed at zero. This maps to OpenCV's CALIB_FIX_K1/K2/K3 flags. When users report "donut artifacts" (peripheral coordinate wraparound), suggest reducing radial_parameters or setting `with_distortion=false` (though the latter is deprecated in current code).

## Common Tasks

### Adding New File Format Support
Extend `io.jl` by:
1. Add detection logic in `load()` function
2. Implement parser that returns a `CalibrationIO` object
3. The existing `Calibration(cio::CalibrationIO)` constructor handles the rest

### Modifying Calibration Pipeline
The pipeline is split across files:
- Corner detection: `_detect_corners()` in detect_fit.jl
- Model fitting: `fit_model()` in detect_fit.jl
- Transformation construction: `obj2img()` and `img2obj()` in buildcalibrations.jl
- Callable interface: `(c::Calibration)(...)` methods in meta.jl

Changes to transformation math typically require updating both forward and inverse pipelines.

### Error Metrics
`calculate_errors()` returns four metrics:
- `reprojection`: Mean pixel error when projecting 3D points to 2D
- `projection`: Mean 3D error when back-projecting 2D points to 3D
- `distance`: Mean error in physical distance between adjacent checkerboard corners
- `inverse`: Mean pixel error after 100 round-trip transformations (2D→3D→2D)

If modifying calibration math, ensure all error metrics remain reasonable (<1 for good calibrations).
