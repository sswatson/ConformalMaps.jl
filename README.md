# ConformalMaps.jl

`ConformalMaps` is a package for approximating the Riemann map from a
simply connected planar domain to a disk. It uses the zipper algorithm as
described in
[Convergence of the zipper algorithm for conformal mapping](http://arxiv.org/abs/math/0605532)
by Don Marshall and Steffen Rohde.

The domain (approximated by a polygon) is specified as an array which lists
the vertices of the domain in counterclockwise order. A conformal
map from `domain` to the unit disk which maps `center` to the origin
is initialized as `ConformalMap(domain,center)`. The keyword argument
`pointspacing=ϵ` inserts equally spaced points along each side of
the polygon so that the spacing between consecutive points is
everywhere less than `ϵ`. Smaller values of `ϵ` give greater accuracy
but require longer to compute. The default value is 1% of the diameter
of the domain. 

```julia
using AsyPlots, ConformalMaps
vertices = [1.0  0.0;
            0.0  1.0;
           -1.0  0.0;
            0.0 -1.0]
f = ConformalMap(vertices,0.0)
```

`f` supports function call syntax: `f(0.1im)`

The inverse of `f` is obtained as `inv(f)` and is of type
`InverseConformalMap`. 

If [`AsyPlots`](https://github.com/sswatson/AsyPlots.jl) is installed,
then `visualize` may be used to display the images of a hyperbolic
tiling of the disk (if called on an `InverseConformalMap`) or grid
lines (if called on a `ConformalMap`).

```julia
g = inv(f)
visualize(g) 
```

![Conformal map](images/square.svg)

`visualize` returns a `ConformalMapVisualization` object, whose
fields `domain` and `range` contain `AsyPlots.Plot2D` objects. 
`combine` returns a single plot with the domain and the range

[![Build Status](https://travis-ci.org/sswatson/ConformalMaps.jl.svg?branch=master)](https://travis-ci.org/sswatson/ConformalMaps.jl)
