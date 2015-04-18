# ConformalMaps.jl

`ConformalMaps` is a package for approximating the Riemann map from a
simply connected planar domain to a disk. It uses the zipper algorithm as
described in
[Convergence of the zipper algorithm for conformal mapping](http://arxiv.org/abs/math/0605532)
by Don Marshall and Steffen Rohde.

The domain (approximated by a polygon) is specified as an array which lists
the vertices of the domain in counterclockwise order. A conformal map is
encoded by an array `ζ` which can be computed with
`initialize_conformal_map`. The keyword argument `resolution=n` inserts `n`
equally spaced points along each side of the polygon. Higher values of `n`
give greater accuracy but require longer to compute. 

```julia
vertices = [1.0  0.0;
	        0.0  1.0;
			-1.0 0.0;
			0.0 -1.0]
ζ = initialize_conformal_map(vertices;resolution=100)
```

Once `ζ` has been computed, a conformal map from the domain to the disk
mapping `center` to the origin can be computed with
`z->conformalmap(ζ,z,center)`. The inverse of this map can be computed with
`z->invconformalmap(ζ,z,center)`.

If [`Graphics2D`](https://github.com/sswatson/Graphics2D.jl) is installed,
then `plotgrid` may be used to display the images of the polar coordinate
level lines under the inverse conformal map. 

```julia
center = 0.0
domain = Complex{Float64}[invconformalmap(ζ,r*exp(im*θ),center) 
                 for r = linspace(0.05,0.99999,20).^(2/3), θ=linspace(0,2π,50)]
showgraphics([plotgrid(domain),Line(closepath(vertices))])
```

![Conformal map](https://github.com/sswatson/ConformalMaps.jl/blob/master/images/conformalmap.png)

[![Build Status](https://travis-ci.org/sswatson/ConformalMaps.jl.svg?branch=master)](https://travis-ci.org/sswatson/ConformalMaps.jl)
