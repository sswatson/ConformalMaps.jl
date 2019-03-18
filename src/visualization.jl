
using AsyPlots
using ConformalMaps

struct ConformalMapVisualization
    domain::AsyPlots.Plot2D
    range::AsyPlots.Plot2D
end

show(io::IO,C::ConformalMapVisualization) = print(io,
                                        "ConformalMapVisualization()")

#-----------------------------------------------------------------------------
# FUNCTIONS FOR DISPLAY VISUALIZING CONFORMAL MAPS
#-----------------------------------------------------------------------------

function closepath(γ::Array{<:Complex,1})
    return closepath(map(AsyPlots.Vec2,γ))
end

function closepath(γ::Array{<:AsyPlots.Vec2})
    if γ[1] == γ[end]
        return γ
    else
        return [γ;[γ[1]]]
    end
end

function makegrid(zvals;ringcolor="blue",raycolor="red",args...)
    lines = Path2D[]
    ringpen = Pen(color=ringcolor)
    for i=1:size(zvals)[1]
        push!(lines,Path2D(hcat(real(zvals[i,:])',
                                imag(zvals[i,:])');pen=ringpen,args...))
    end
    raypen = Pen(color=raycolor)
    for j=1:size(zvals)[2]
        push!(lines,Path2D(hcat(real(zvals[:,j]),
                                imag(zvals[:,j]));pen=raypen,args...))
    end
    return lines
end

function hyperbolictiling(f::Union{Function,InverseConformalMap};
                          rings::Integer=9,
                          rays::Integer=16,
                          rotation::Real=0.0,
                          innerradius::Real=1.0/3.0,
                          lwfunction=(i->2^(-i/2)), 
                          ringcolor="blue",
                          raycolor="red",
                          spline=true)
    
   points = [[f((1-(1-innerradius)/2^(k-1))*cos(θ+rotation) +
                  im*(1-(1-innerradius)/2^(k-1))*sin(θ+rotation))
                       for θ=range(0, stop=2π, length=1+rays*2^(k-1))] for k=1:rings]
   return Plot([
            (isa(f,InverseConformalMap) ? domain(f) :
             [Point2D(f(0.0),linewidth=2),
              Polygon2D([f(cis(θ)) for θ=range(0, stop=2π, length=250)];linewidth=0.3)]);
             [Path2D([points[i][k],points[i+1][2*k-1]];
                     color=raycolor,linewidth=lwfunction(i))
                          for i=1:length(points)-1 for k=1:length(points[i])-1]
             [Polygon2D(points[i][1:end-1];
                    color=ringcolor,linewidth=lwfunction(i),spline=spline)
                    for i=1:length(points)]])
end

function domain(f::Union{ConformalMap,InverseConformalMap})
    return [Point2D(f.center,linewidth=2),Path2D(closepath(f.domain),linewidth=0.3)]
end

"""
    intersectQ(p,q,r,s)

Determine whether the segment [p,q] intersects
the segment [r,s].

The four points can be 2-tuples of `Real`s, or
`Complex`es.

    intersectQ(p,q,γ)

Determine whether the segment from p to q intersects
the polygonal path γ
"""
function intersectQ(p::Tuple{<:Real,<:Real},
                    q::Tuple{<:Real,<:Real},
                    r::Tuple{<:Real,<:Real},
                    s::Tuple{<:Real,<:Real})
    a = -((p[2]*(r[1] - s[1]) + r[2]*s[1] - r[1]*s[2] +
          p[1]*(-r[2] + s[2]))/((-p[2] + q[2])*(r[1] - s[1]) + (p[1] -
                 q[1])*(r[2] - s[2])))
    b = (-q[2]*r[1] + p[2]*(-q[1] + r[1]) + p[1]*(q[2] - r[2]) +
         q[1]*r[2])/((p[2] - q[2])*(r[1] - s[1]) - (p[1] - q[1])*(r[2] - s[2]))
    return 0 <= a <= 1 && 0 <= b <= 1
end

function intersectQ(p::Complex,q::Complex,r::Complex,s::Complex)
    # Determines whether the segment from p to q intersects the
    # segment from r to s
    p,q,r,s = map(reim,(p,q,r,s))
    return intersectQ(p,q,r,s)
end

function intersectQ(p::Complex,q::Complex,gamma::Array{<:Complex,1})
    for i=1:length(gamma)-1
        if intersectQ(p,q,gamma[i],gamma[i+1])
            return true
        end
    end
    return false
end

function intersectQ(p::AsyPlots.Vec2,
                    q::AsyPlots.Vec2,
                    γ::Array{<:AsyPlots.Vec2,1})
    intersectQ(map(complex,(p,q))...,map(complex,γ))
end

function makegrid(boundary::Array{<:AsyPlots.Vec2,1},n::Integer)
    grid  = Tuple[]
    xvals = [P.x for P in boundary]
    yvals = [P.y for P in boundary]
    ϵ = max(maximum(xvals)-minimum(xvals),maximum(yvals)-minimum(yvals))/(n-1)
    m = length(minimum(xvals):ϵ:maximum(xvals))
    n = length(minimum(yvals):ϵ:maximum(yvals))
    totalgrid = [(x,y) for
        x = range(minimum(xvals), stop=maximum(xvals), length=m),
        y = range(minimum(yvals), stop=maximum(yvals), length=n)]
    pointsinside = [iswellinside(AsyPlots.Vec2(totalgrid[i,j]),
                                boundary;epsilon=1e-3) ? 1 : 0 for i=1:m,j=1:n]
    lines = Tuple{Tuple{Int64,Int64},Tuple{Int64,Int64}}[]
    for i=1:m
        for j=1:n
            for v in [(i+1,j),(i-1,j),(i,j+1),(i,j-1)]
                if 1 <= v[1] <= m && 1 <= v[2] <= n &&
                    pointsinside[i,j] == 1 && pointsinside[v...] == 1 &&
                       ~intersectQ(AsyPlots.Vec2(totalgrid[i,j]),
                                   AsyPlots.Vec2(totalgrid[v...]),
                                   boundary)
                    push!(lines,((i,j),v))
                end
            end
        end
    end
    return totalgrid, pointsinside, lines
end

makegrid(boundary::Array{T,1},n::Integer) where {T<:Complex} =
    makegrid(hcat(reim(boundary)...),n)

function grid(domain,totalgrid,pointsinside,lines,center)
    return AsyPlots.Plot([[AsyPlots.Point(real(center),imag(center)),
    AsyPlots.Path(closepath(domain))];
    [AsyPlots.Path([totalgrid[line[1]...],
                    totalgrid[line[2]...]],linewidth=0.2) for line in lines]])
end

function gridimage(f,totalgrid,pointsinside,lines,center)
    return Plot([[AsyPlots.Point(0.0,0.0),AsyPlots.Circle((0.0,0.0),1,linewidth=0.3)];
     [AsyPlots.Path([f(totalgrid[line[1]...]),
                     f(totalgrid[line[2]...])],linewidth=0.1) for line in lines]])
end

function visualize(CM::ConformalMap,n::Integer=25)
    totalgrid,pointsinside,lines = makegrid(closepath(CM.domain),n)
    return ConformalMapVisualization(
        grid(closepath(CM.domain),totalgrid,pointsinside,lines,CM.center),
        gridimage(p->CM(p[1] + im*p[2]),totalgrid,pointsinside,lines,CM.center)
    )
end

function visualize(ICM::InverseConformalMap;kwargs...)
    return ConformalMapVisualization(
        hyperbolictiling(z->z;kwargs...),hyperbolictiling(ICM;kwargs...)
    )
end

function combine(V::ConformalMapVisualization;kwargs...)
    domainbb = AsyPlots.boundingbox(V.domain)
    rangebb = AsyPlots.boundingbox(V.range)
    totalwidth = domainbb.xmax - domainbb.xmin + rangebb.xmax - rangebb.xmin 
    xshift = domainbb.xmax - rangebb.xmin + 0.2*totalwidth
    domainycenter = mean([domainbb.ymax,domainbb.ymin])
    rangeycenter = mean([rangebb.ymax,rangebb.ymin])
    yshift = domainycenter - rangeycenter 
    fullpicture = V.domain + Shift(xshift,yshift)*V.range
    bb = AsyPlots.boundingbox(fullpicture)
    midheight = mean([bb.ymax,bb.ymin]) 
    arrow = Path([domainbb.xmax+0.04*totalwidth midheight;
                  domainbb.xmax+0.16*totalwidth midheight];arrow=AsyPlots.Arrow(4))
    fullpicture += Plot(arrow)
    Plot(fullpicture.elements;kwargs...) 
end

function Base.show(io::IO,mime::MIME"text/plain",V::ConformalMapVisualization)
    Base.show(io,mime,combine(V))
end

function Base.show(io::IO, mime::MIME"image/svg+xml", V::ConformalMapVisualization)
    Base.show(io,mime,combine(V))
end

function Base.show(io::IO, mime::MIME"image/png", V::ConformalMapVisualization)
    Base.show(io,mime,combine(V))
end

#-----------------------------------------------------------------------------
