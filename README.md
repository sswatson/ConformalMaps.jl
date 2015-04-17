# ConformalMaps.jl

[![Build Status](https://travis-ci.org/sswatson/ConformalMaps.jl.svg?branch=master)](https://travis-ci.org/sswatson/ConformalMaps.jl)

`ConformalMaps` is a package for approximating the Riemann map from a
simply connected planar domain to a disk. It uses the zipper algorithm as
described in
[Convergence of the zipper algorithm for conformal mapping](http://arxiv.org/abs/math/0605532)
by Don Marshall and Steffen Rohde.

The domain (approximated by a polygon) is specified as an array which lists
the vertices of the domain in counterclockwise order. A conformal map is
encoded by an array `zeta` which can be computed with
`initialize_conformal_map`. The keyword argument `resolution=n` inserts `n`
equally spaced points along each side of the polygon. Higher values of `n`
give greater accuracy but require longer to compute. 

```julia
zeta = initialize_conformal_map([1.0  0.0;
	                             0.0  1.0;
							    -1.0 0.0;
								 0.0 -1.0];resolution=100)
```

Once `zeta` has been computed, a conformal map from the domain to the disk
mapping `center` to the origin can be computed with
`z->conformalmap(zeta,z,center)`. The inverse of this map can be computed with
`z->invconformalmap(zeta,z,center)`.

If [`Graphics2D`](https://github.com/sswatson/Graphics2D.jl) is installed,
then `plotgrid` may be used to display the images of the polar coordinate
level lines under the inverse conformal map. 

```julia
domain = Complex{Float64}[invconformalmap(zeta,r*exp(im*theta),center) 
                 for r=linspace(0,0.99999,100),theta=linspace(0,2*pi,100)]
showgraphics([plotgrid(domain),Line(closepath(mypoints))])
```

![Conformal map](https://github.com/sswatson/ConformalMaps.jl/images/square.png)

