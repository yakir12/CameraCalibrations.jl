using ImageDraw
dir = joinpath(@__DIR__(), "example")
files = filter(file -> last(splitext(file)) == ".png", readdir(dir, join = true))
n = length(files)
n_corners = (5, 8)
d = Dict()
for file in files
    res = CameraCalibrations._detect_corners(file, n_corners)
    if !ismissing(res) 
        _, xys = res
        img = RGB.(FileIO.load(file))
        for xy in xys
            draw!(img, Cross(Point(round.(Int, reverse(xy))...), 5), RGB(1,0,0))
        end
        FileIO.save(basename(file), img)
        d[basename(file)] = vec(xys)
    end
end
open("example/corners.json", "w") do io
    print(io, JSON3.write(d))
end

