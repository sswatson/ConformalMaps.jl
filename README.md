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
`resolution=n` inserts `n` equally spaced points along each side of
the polygon. Higher values of `n` give greater accuracy but require
longer to compute. 

```julia
using Graphics2D
using ConformalMaps
vertices = [1.0  0.0;
            0.0  1.0;
           -1.0  0.0;
            0.0 -1.0]
f = ConformalMap(vertices,0.0;resolution=100)
```

`f` supports function call syntax: `f(0.1im)`

The inverse of `f` is obtained as `inv(g)` and is of type
`InverseConformalMap`. 

If [`Graphics2D`](https://github.com/sswatson/Graphics2D.jl) is installed,
then `visualize` may be used to display the images of certain polar
coordinate level lines under the inverse conformal map (if called on
an `InverseConformalMap`) or grid lines (if called on a `ConformalMap`). 

```julia
g = inv(f)
showgraphics([visualize(g;rays=24,rings=12,innerradius=0.2)[2];domain(f)])
```

![Conformal map](https://github.com/sswatson/ConformalMaps.jl/blob/master/images/square.png)

[![Build Status](https://travis-ci.org/sswatson/ConformalMaps.jl.svg?branch=master)](https://travis-ci.org/sswatson/ConformalMaps.jl)
