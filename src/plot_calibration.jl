function get_axes(ratio, checker_size, n_corners, sz)
    w = round.(Int, ratio .* checker_size .* (n_corners .- 1))
    min_col, min_row = round.(Int, (w .- sz) ./ 2)
    max_col, max_row = (min_col, min_row) .+ sz .- 1
    return min_col:max_col, min_row:max_row
end

function get_ratio(imgpoints, checker_size)
    l = mean(1:2) do dims
        mean(norm, diff(imgpoints; dims))
    end
    l/checker_size
end

function image_transformations(c, extrinsic_index, imgpointss, checker_size, n_corners, sz)
    ratio = get_ratio(imgpointss[extrinsic_index], checker_size)
    s = LinearMap(SDiagonal{2}(ratio*I))
    tform = c.real2image[extrinsic_index] ∘ Base.Fix2(push, 0) ∘ inv(s)
    itform = s ∘ pop ∘ c.image2real[extrinsic_index]
    axs = get_axes(ratio, checker_size, n_corners, sz)
    return (; tform, itform, axs)
end

function draw_crosses!(img, imgpoint, n1, color)
    ij = vec([round.(Int, rc) for rc in imgpoint])
    radius = round(Int, norm(ij[1] - ij[n1])/n1/5) # radius relative to checkerboard's aparent size
    draw!(img, Cross.(Point.(CartesianIndex.(Tuple.(ij))), radius), color)
end

plot(::Nothing, c::Calibration, imgpointss, n_corners, checker_size, sz) = nothing

"""
    plot
Save to locally created directory `debug` all the calibration's images after rectification. Blue crosses indicate detetcted corners, and red crosses indicate their rectificated real-world locations. Ideally, both crosses' centers are right on top of each other.
"""
function plot(plot_folder, c::Calibration, imgpointss, n_corners, checker_size, sz)
    dir = mkpath(plot_folder)
    for (file, imgpoints) in zip(c.files, imgpointss)
        img = RGB.(FileIO.load(file))
        # img = RGB.(permutedims(FileIO.load(file)))
        draw_crosses!(img, imgpoints, n_corners[1], colorant"red")
        tform, itform, axs = image_transformations(c, findfirst(==(file), c.files), imgpointss, checker_size, n_corners, sz)
        imgw = warp(img, tform, axs)
        draw_crosses!(imgw, itform.(imgpoints), n_corners[1], colorant"blue")
        FileIO.save(joinpath(dir, basename(file)), parent(imgw))
    end
    return nothing
end
