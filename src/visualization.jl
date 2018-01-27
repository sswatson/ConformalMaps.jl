
using Graphics2D
using ConformalMaps

immutable ConformalMapVisualization
    domain::Array{<:Graphics2D.GraphicElement,1}
    range::Array{<:Graphics2D.GraphicElement,1}
end

show(io::IO,C::ConformalMapVisualization) = print(io,"ConformalMapVisualization()")

#-----------------------------------------------------------------------------
# FUNCTIONS FOR DISPLAY VISUALIZING CONFORMAL MAPS 
#-----------------------------------------------------------------------------

function closepath{T<:AbstractFloat}(γ::Array{Complex{T},1})
    return closepath(hcat([real(a) for a in γ], [imag(a) for a in γ]))
end

function closepath{T<:Real}(γ::Array{Complex{T},1})
    return closepath(float(γ))
end

function closepath{T<:AbstractFloat}(γ::Array{T,2})
    if γ[end,:] == γ[1,:]
        return γ
    else
        return vcat(γ,γ[1,:]')
    end
end

function closepath{T<:Real}(γ::Array{T,2})
    return closepath(float(γ))
end

function makegrid(zvals;color1="blue",color2="red",args...)
    lines = Graphics2D.Line[]
    for i=1:size(zvals)[1]
        push!(lines,Graphics2D.Line(hcat(real(zvals[i,:])',
                                         imag(zvals[i,:])');color=color1,args...))
    end
    for j=1:size(zvals)[2]
        push!(lines,Graphics2D.Line(hcat(real(zvals[:,j]),
                                         imag(zvals[:,j]));color=color2,args...))
    end
    return lines
end

function hyperbolictiling(f::Function;
                          rings::Integer=9,
                          rays::Integer=16,
                          rotation::Real=0.0,
                          innerradius::Real=1.0/3.0,
                          ringcolor="blue",
                          raycolor="red")
    points = Array{Complex64,1}[[f((1-(1-innerradius)/2^(k-1))*cos(θ+rotation) + 
                                   im*(1-(1-innerradius)/2^(k-1))*sin(θ+rotation)) 
                                 for θ=linspace(0,2π,1+rays*2^(k-1))] for k=1:rings]
    return [vcat([[Graphics2D.Line([points[i][k],points[i][k+1]];
                color=ringcolor,linewidth=1-i/(rings+4)) 
            for k=1:length(points[i])-1] 
            for i=1:length(points)]...);
            vcat([[Graphics2D.Line([points[i][k],points[i+1][2*k-1]];
                color=raycolor,linewidth=1-i/(rings+4)) 
            for k=1:length(points[i])-1] 
            for i=1:length(points)-1]...)]
end

function distance(p::Tuple{<:Real,<:Real},q::Array{<:Real,2},r::Array{<:Real,2})
    if r[1] == q[1]
        if min(r[2],q[2]) < p[2] < max(r[2],q[2])
            return abs(p[1]-r[1])
        else
            return min(hypot(p[1]-q[1],p[2]-q[2]),hypot(p[1]-r[1],p[2]-r[2]))
        end
    else
        m = (r[2]-q[2])/(r[1]-q[1])
        x = (q[1]*m^2 - q[2]*m + p[2]*m + p[1])/(m^2+1)
        y = (p[2]*m^2 - q[1]*m + p[1]*m + q[2])/(m^2+1)
        if min(r[1],q[1]) < x < max(r[1],q[1])
            return hypot(p[1]-x,p[2]-y)
        else
            return min(hypot(p[1]-q[1],p[2]-q[2]),hypot(p[1]-r[1],p[2]-r[2]))
        end
    end
end

function distance(p::Tuple{<:Real,<:Real},gamma::Array{<:Real,2})
    m = hypot(p[1]-gamma[1,1],p[2]-gamma[1,2])
    for i=1:size(gamma)[1]-1
        d = distance(p,gamma[i,:],gamma[i+1,:])
        if d < m
            m = d
        end
    end
    return m
end

function inside(p::Tuple{<:Real,<:Real},gamma::Array{<:Real,2},tol::Float64=1e-8)
    if gamma[1,:] != gamma[length(gamma[:,1]),:] 
        return false
    end
    if distance(p,gamma) < 1e-3
        return false
    end
    cntr = 0; m = sqrt(2); # the slope is an arbitrary irrational number
    for i=1:length(gamma[:,1])-1
        (x1,y1,x2,y2) = (gamma[i,1],gamma[i,2],gamma[i+1,1],gamma[i+1,2])
        if ((y2 - p[2] - m*(x2-p[1]))*(y1 - p[2] - m*(x1-p[1])) < 0) 
            if (m*p[1]*x1 - p[2]*x1 - m*p[1]*x2 + p[2]*x2 - x2*y1 + x1*y2)/(m*x1 - m*x2 - y1 + y2)  - p[1] > 0
                cntr += 1
            end
        end
    end
    return isodd(cntr)    
end

inside(p::Tuple{<:Real,<:Real},gamma::Array{<:Complex,1},tol::Float64=1e-8) = 
    inside(p,hcat(reim(gamma)...),tol)

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

function intersectQ(p::Tuple{<:Real,<:Real},
                   q::Tuple{<:Real,<:Real},
                   gamma::Array{<:Real,2})
    for i=1:size(gamma)[1]-1
        if intersectQ(p,q,(gamma[i,1],gamma[i,2]),(gamma[i+1,1],gamma[i+1,2]))
            return true
        end
    end
    return false
end

function domain(f::Union{ConformalMap,InverseConformalMap})
    return Graphics2D.GraphicElement[Graphics2D.Point(f.center),
                                     Graphics2D.Line(closepath(f.domain))]
end

function makegrid(boundary::Array{T,2} where T<:Real,n::Integer)
    grid = Tuple{T,T}[]
    xvals = boundary[:,1]
    yvals = boundary[:,2]
    ϵ = max(maximum(xvals)-minimum(xvals),maximum(yvals)-minimum(yvals))/(n-1)
    m = length(minimum(xvals):ϵ:maximum(xvals))
    n = length(minimum(yvals):ϵ:maximum(yvals))
    totalgrid = [(x,y) for 
        x = linspace(minimum(xvals),maximum(xvals),m),
        y = linspace(minimum(yvals),maximum(yvals),n)]
    pointsinside = zeros(Int64,m,n)
    for i=1:m
        for j=1:n
            if inside(totalgrid[i,j],boundary)
                pointsinside[i,j] = 1
            end
        end
    end
    lines = Tuple{Tuple{Int64,Int64},Tuple{Int64,Int64}}[]
    for i=1:m
        for j=1:n 
            for v in [(i+1,j),(i-1,j),(i,j+1),(i,j-1)]
                if 1 <= v[1] <= m && 1 <= v[2] <= n && 
                    pointsinside[i,j] == 1 && pointsinside[v...] == 1 && 
                       ~intersectQ(totalgrid[i,j],totalgrid[v...],boundary)
                    push!(lines,((i,j),v))
                end
            end
        end
    end
    return totalgrid, pointsinside, lines
end

makegrid{T<:Complex}(boundary::Array{T,1},n::Integer) = 
    makegrid(hcat(reim(boundary)...),n)

function grid(domain,totalgrid,pointsinside,lines,center)
    return [[Graphics2D.Point(real(center),imag(center)),
    Graphics2D.Line(closepath(domain))];
    [Graphics2D.Line([totalgrid[line[1]...],
    totalgrid[line[2]...]],linewidth=0.2) for line in lines]] 
end

function gridimage(f,totalgrid,pointsinside,lines,center)
    return [[Graphics2D.Point(0.0,0.0),Graphics2D.Circle([0.0,0.0],1)];
     [Graphics2D.Line([f(totalgrid[line[1]...]),
                       f(totalgrid[line[2]...])],linewidth=0.1) for line in lines]] 
end

function visualize(CM::ConformalMap,n::Integer=40)
    totalgrid,pointsinside,lines = makegrid(closepath(CM.domain),n)
    return ConformalMapVisualization(
        grid(closepath(CM.domain),totalgrid,pointsinside,lines,CM.center), 
        gridimage(p->CM(p[1] + im*p[2]),totalgrid,pointsinside,lines,CM.center)
    )
end

function visualize(ICM::InverseConformalMap;kwargs...)
    return ConformalMapVisualization(
        hyperbolictiling(z->z;kwargs...),hyperbolictiling(z->ICM(z);kwargs...)
    )   
end

#-----------------------------------------------------------------------------
