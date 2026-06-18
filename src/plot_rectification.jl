"""
    reference_affine(c)
The affine (planar) part of `rectification(c)`: the `Translation` to the `center` origin and, when a `north` is set,
the rotation aligning `center → north` to `+y`. Equivalent to `rectification(c)` with the leading `pop ∘
image2real[extrinsic_index]` factored out, so that `rectification(c) == reference_affine(c) ∘ pop ∘
image2real[extrinsic_index]`. Returned as an invertible `CoordinateTransformations` map.
"""
function reference_affine(c::Calibration)
    f = pop ∘ c.image2real[c.extrinsic_index]
    cc = f(c.center)
    A = Translation(-cc)
    isnothing(c.north) && return A
    p = f(c.north) - cc
    θ = π/2 - atan(p[2], p[1])
    return LinearMap(Angle2d(θ)) ∘ A
end

plot_rectification(::Nothing, c::Calibration, imgpointss, n_corners, checker_size) = nothing

"""
    plot_rectification
Save `rectified_extrinsic.png` to `plot_folder`: the chosen extrinsic image (`c.files[c.extrinsic_index]`) warped
into the reference frame defined by `c.center`/`c.north`. Use it to confirm the `center` lands at the image center
and `north` points up. The image is rendered `+y`-up (north up, `+x` to the right) so it reads like the original
photo, without mirroring. Requires `c.center` and `c.extrinsic_index` to be set.
"""
function plot_rectification(plot_folder, c::Calibration, imgpointss, n_corners, checker_size)
    ei = c.extrinsic_index
    ratio = get_ratio(imgpointss[ei], checker_size)
    # display map: reference-frame xy → output pixel (row, col), `+y`-up and scaled to ~`ratio` pixels per checker.
    # Angle2d(π/2) sends (x, y) → (-y, x), so row = -ratio*y (north points up) and col = ratio*x (+x to the right);
    # a proper rotation (det > 0), hence no mirror.
    D = LinearMap(Angle2d(π/2)) ∘ LinearMap(SDiagonal{2}(ratio*I))
    A = reference_affine(c)
    M = D ∘ A                                   # world xy → output pixel
    P = rectification(c)                        # source pixel → reference-frame xy

    img = RGB.(FileIO.load(c.files[ei]))
    nrow, ncol = size(img)
    # output extents: bounding box of the four source-image corners mapped into the output frame.
    corners = (RowCol(1, 1), RowCol(1, ncol), RowCol(nrow, 1), RowCol(nrow, ncol))
    ds = map(D ∘ P, corners)
    rows = [d[1] for d in ds]
    cols = [d[2] for d in ds]
    rowrange = floor(Int, minimum(rows)):ceil(Int, maximum(rows))
    colrange = floor(Int, minimum(cols)):ceil(Int, maximum(cols))

    tform = c.real2image[ei] ∘ Base.Fix2(push, 0) ∘ inv(M)
    imgw = warp(img, tform, (rowrange, colrange))
    FileIO.save(joinpath(mkpath(plot_folder), "rectified_extrinsic.png"), parent(imgw))
    return nothing
end
