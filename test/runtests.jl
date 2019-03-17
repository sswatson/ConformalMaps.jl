using ConformalMaps
using Test
using Requires 

vertices = [1.0  0.0;
            0.0  1.0;
            -1.0  0.0;
            0.0 -1.0]
f = ConformalMap(vertices,0.0)
@test f(0.0) == 0

vertices = [1, im, -1, -im] 
f = ConformalMap(vertices,0.0im)
g = inv(f)
z = 0.1 + 0.1im; 
@test g(f(z)) ≈ z 

vertices = Complex{BigFloat}[1, im, -1, -im] 
f = ConformalMap(vertices,0.0im)
g = inv(f)
z = 0.1 + 0.1im; 
@test g(f(z)) ≈ z 

@require AsyPlots begin
    slitdomain = [0.0, 0.495, 0.5 + 0.25*im, 0.505, 1.0, 1.0 + 1.0*im, im]
    f = ConformalMap(slitdomain,0.5+0.5im;resolution=40)
    g = inv(f)
    visualize(g)
end
