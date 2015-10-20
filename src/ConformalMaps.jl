module ConformalMaps

# This package provides a ConfomalMap type which numerically approximates the 
# conformal map from a Jordan domain specified by a list of points on the
# boundary to a disk or half-plane. 

# The algorithm used is based on the paper Convergence of the Zipper algorithm 
# for conformal mapping by Don Marshall and Steffen Rohde. 
# http://arxiv.org/abs/math/0605532

#-----------------------------------------------------------------------------
# BASIC USAGE 
#-----------------------------------------------------------------------------

## Computing the list of points ζ can take a long time for a complicated domain
# f = ConformalMap(Float64[1 0; 0 1; -1 0; 0 -1];resolution=100)

## ConformalMap types support function call notation
# f(0.6+0.7im)

## The inverse of a ConformalMap is an InverseConformalMap
# g = inv(f)

## visualize(g::InverseConformalMap) outputs two graphics, the second
## of which shows the images under g of a hyperbolic tiling of the disk. 
# showgraphics(visualize(g;rays=16,rings=12)[2])

#-----------------------------------------------------------------------------


import Graphics2D,
       Base.call,
       Base.show,
       Base.inv,
       Base.intersect

export ConformalMap,
       closepath,
       D_to_H,
       H_to_D,
       hyperbolictiling,
       visualize,
       domain

#-----------------------------------------------------------------------------
# CONFORMAL MAP IMMUTABLE TYPE
#-----------------------------------------------------------------------------

immutable ConformalMap{T<:Real} # A map from the domain to the disk
    domain::Array{Complex{T},1} # Points tracing out the boundary
    data::Array{Complex{T},1} # An array (ζ in the paper) encoding the conformal map
    center::Complex{T} # The point in the domain which maps to the origin
end

immutable InverseConformalMap{T<:Real} # A map from the disk to the domain
    domain::Array{Complex{T},1} 
    data::Array{Complex{T},1} 
    center::Complex{T} 
end

function ConformalMap{T<:Real}(domain::Array{Complex{T},1},
                               center::Complex{T};
                               resolution::Integer=1)

    densepoints = densify(domain,resolution)

    ζ = Complex{typeof(real(domain[1]))}[]

    for k=3:length(densepoints)
        push!(ζ,fcompose(ζ,phi2(densepoints[k],
                                densepoints[1],densepoints[2])))
    end

    push!(ζ,fcompose(ζ,phi2(densepoints[1],densepoints[1],densepoints[2])))

    push!(ζ,densepoints[1])
    push!(ζ,densepoints[2])

    push!(ζ,philast(fcompose(ζ[1:end-3],phi2(center,ζ[end-1],ζ[end])),ζ[end-2]))

    return ConformalMap(densepoints,ζ,center)

end

function ConformalMap{T<:Real}(domain::Array{Complex{T},1},
                       center::T;
                       kwargs...)
    return ConformalMap(domain,complex(center);kwargs...)
end

function ConformalMap{T<:Real}(domain::Array{T,2},
                               center::T;
                               kwargs...)
    return ConformalMap(domain[:,1] + im*domain[:,2],center;kwargs...)
end

function InverseConformalMap{T<:Real}(domain::Array{Complex{T},1},
                       center::Complex{T};kwargs...)
    return inv(ConformalMap(domain,center;kwargs...))
end

function inv(CM::ConformalMap)
    return InverseConformalMap(CM.domain,CM.data,CM.center)
end

function inv(ICM::InverseConformalMap)
    return ConformalMap(ICM.domain,ICM.data,ICM.center)
end

function Base.call(CM::ConformalMap,z::Union{Real,Complex})
    ζ = CM.data
    return H_to_D((philast(fcompose(ζ[1:end-4],
                                    phi2(z,ζ[end-2],ζ[end-1])),ζ[end-3])),ζ[end])
end

function Base.call(ICM::InverseConformalMap,z::Union{Real,Complex})
    ζ = ICM.data
    return phi2inv(finvcompose(ζ[1:end-4],
                               philastinv(D_to_H(z,ζ[end]),ζ[end-3])),ζ[end-2],ζ[end-1])
end

function shortformat(x::AbstractFloat)
    if x == 0
        return "0"
    elseif abs(x) > 0.1
        return @sprintf("%.2f",x)
    elseif abs(x) > 0.01
        return @sprintf("%.3f",x)
    elseif abs(x) > 0.001
        return @sprintf("%.4f",x)
    else
        return @sprintf("%.2e",x)
    end
end

function shortformat(z::Complex)
    if imag(z) == 0
        return shortformat(real(z))
    elseif real(z) == 0
        return string(shortformat(imag(z)),"im")
    else
        return string(shortformat(real(z)), "+", shortformat(imag(z)),"im")
    end
end

function show(io::IO,CM::Union{ConformalMap,InverseConformalMap})
    if typeof(CM) <: ConformalMap
        print(io,"ConformalMap{")
    else
        print(io,"InverseConformalMap{")
    end
    print(io,typeof(CM).parameters[1])
    print(io,"}([")
    print(io,shortformat(CM.domain[1]))
    print(io,",")
    print(io,shortformat(CM.domain[2]))
    print(io,",...,")
    print(io,shortformat(CM.domain[end-1]))
    print(io,",")
    print(io,shortformat(CM.domain[end]))
    print(io,"],center=")
    print(io,shortformat(CM.center))
    print(io,")")
end

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# COMPLEX NUMBER CONVERSION
#-----------------------------------------------------------------------------
# for conveniently splitting complex numbers and arrays of them into 
# real and imaginary parts 
realify{T<:Complex}(A::Array{T,1}) = hcat(real(A),imag(A))
realify{T<:AbstractFloat}(z::Complex{T}) = (real(z),imag(z))
realify(x::Real) = (x,0.0)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# EXPLICIT CONFORMAL MAPS
#-----------------------------------------------------------------------------
# f and finv are the basic slit-domain conformal maps composed to 
# construct the global conformal map 
function f{T<:AbstractFloat}(z::Union{T,Complex{T}},
                             a::Union{T,Complex{T}};
                             tol::Float64=1e-12)
    if isnan(z)
        return -im*sign(imag(a))*abs(a)^2*sqrt(1/real(a)^2+1/imag(a)^2)
    elseif abs(imag(a)) < tol
        return sqrt(z^2 - a^2)
    else
        return sqrt((z/(1-z/(im*abs(a)^2/imag(a))))^2-abs(a)^4/real(a)^2)
    end
end

function finv{T<:AbstractFloat}(w::Union{T,Complex{T}},
                                a::Union{T,Complex{T}};
                                tol::Float64=1e-12)
    if abs(imag(a)) < tol
        return sqrt(w^2+a^2)
    elseif isnan(w)
        return im*abs(a)^2/imag(a)
    else
        return abs(a)^2/imag(a)*sqrt(w^2+(abs(a)^2/real(a))^2) / 
               (abs(a)^2/imag(a)-im*sqrt(w^2+(abs(a)^2/real(a))^2))
    end
end

function fcompose{T<:AbstractFloat}(A::Array{Complex{T},1},
                                    z::Union{AbstractFloat,Complex})
    for i=1:length(A)
        z = f(z,A[i])
    end
    return z
end

function finvcompose{T<:AbstractFloat}(A::Array{Complex{T},1},
                                       w::Union{T,Complex{T}})
    for i=length(A):-1:1
        w = finv(w,A[i])
    end
    return w
end

# The first and last maps (phi2 and philast) in Figure 3 in [MR06] have special forms

function phi2(z::Union{Complex,AbstractFloat},
              point1::Union{Complex,AbstractFloat},
              point2::Union{Complex,AbstractFloat})
    if z == point1
        return convert(typeof(z),NaN)
    else
        return sqrt((z-point2)/(z-point1))
    end
end

phi2inv(w::Union{Complex,AbstractFloat},
        point1::Union{Complex,AbstractFloat},
        point2::Union{Complex,AbstractFloat}) = (w^2*point1-point2)/(w^2-1.0)
philast(z::Union{Complex,AbstractFloat},
        ζlast::Union{Complex,AbstractFloat}) = -(z/(1-z/(ζlast)))^2
philastinv(w::Union{Complex,AbstractFloat},
           ζlast::Union{Complex,AbstractFloat}) = (w*ζlast+im*sqrt(w)*ζlast^2)/(w+ζlast^2)

# Explicit conformal maps between the upper half H and the disk D
H_to_D(z::Union{Complex,AbstractFloat},
       a::Union{Complex,AbstractFloat}) = (z-a)/(z-conj(a))

D_to_H(w::Union{Complex,AbstractFloat},
       a::Union{Complex,AbstractFloat}) = (w*conj(a)-a)/(w-1);
#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------
# GEODESIC ALGORITHM CONFORMAL MAPS
#-----------------------------------------------------------------------------

# This function adds additional points along the segments defining 
# the boundary of the domain
function densify{T<:AbstractFloat}(A::Array{T,2},n::Integer)
    if n == 1
        return A
    end
    densearray = zeros(typeof(A[1,1]),n*size(A)[1],2)
    for j=1:size(A)[1]-1
        for k=1:n
            densearray[(j-1)*n+k,:] = (n+1-k)/(n)*A[j,:] + (k-1)/(n)*A[j+1,:]
        end
    end
    for k=1:n
        densearray[(size(A)[1]-1)*n+k,:] = (n+1-k)/(n)*A[size(A)[1],:] + (k-1)/(n)*A[1,:]
    end
    return densearray
end

function densify{T<:AbstractFloat}(A::Array{Complex{T},1},n::Integer)
    if n == 1
        return A
    end
    densearray = zeros(typeof(A[1]),n*length(A))
    for j = 1:length(A)-1
        for k = 1:n
            densearray[(j-1)*n+k] = (n+1-k)/n * A[j] + (k-1)/n * A[j+1]
        end
    end
    for k=1:n
        densearray[(length(A)-1)*n+k] = (n+1-k)/(n)*A[length(A)] + (k-1)/(n)*A[1]
    end
    return densearray
end

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# ZIPPER ALGORITHM FUNCTIONS
#-----------------------------------------------------------------------------

function newton(f::Function,
                fprime::Function,
                x0::Union{Real,Complex};
                tol::Float64=1e-15,
                maxiter::Int64=100,
                UHP::Bool=false)
    function reflect(xnew::Union{Real,Complex})
        if UHP && imag(xnew) < 0.0
            xnew = real(xnew) + 0.0*im
        end
        return xnew
    end
    xold = x0
    xnew = xold - f(xold)/fprime(xold)
    xnew = reflect(xnew)
    cntr = 1
    while abs(f(xnew)) > tol
        xold = xnew
        xnew = xold - f(xold)/fprime(xold)
        xnew = reflect(xnew)
        cntr += 1
        if cntr > maxiter
            return xnew # error("Maximum iterations reached in Newton's method: $xnew") 
        end
    end
    return xnew
end

function fzip(d::Complex,z::Union{Real,Complex})
    p = angle(d)/π
    return abs(d)/(p^p*(1-p)^(1-p)) * (z-p)^(p)*(z-(p-1))^(1-p)
end

function fzip_pr(d::Complex,z::Union{Real,Complex})
    p = angle(d)/π
    return abs(d)/(p^p*(1-p)^(1-p)) * 
    p*(z-p)^(p-1)*(1-p+z)^(1-p) + (1-p)*(z-p)^p*(1-p+z)^(-p)
end

function fzipinv(a::Complex,w::Union{Real,Complex};args...)
    p = angle(a)/π
    d = exp(im*p*π)*p^p*(1-p)^(1-p)
    v = complex(w / (abs(a)/(p^p*(1-p)^(1-p))))
    if abs(v) > 9/8 * abs(d)
        return newton(z->fzip(d,z)-v,z->fzip_pr(d,z), v + 2p - 1 + p*(1-p)/(2*v) + 
                                (1 - 2p)*p*(1-p)/(3*v*v);args...)
    elseif abs(v - d) < 0.25 * imag(d)
        k(z) = im*sqrt((z-d)*exp(-im*p*π))
        kpr(z) = im*exp(-im*π*p)/sqrt((z-d)*exp(-im*π*p))
        return newton(z->k(fzip(d,z))-k(v),z->kpr(fzip(d,z))*fzip_pr(d,z),k(v);args...)
    elseif angle(v) < π*p
        kofzip(z) = (z-p)*((z-(p-1)))^(1/p-1)
        kof_prime(z) = (z-p)*(1/p-1)*(z-(p-1))^(1/p-2) + (z-(p-1))^(1/p-1)
        u = v^(1/p)
        return newton(z->kofzip(z)-u,kof_prime,p+u;args...)
    else 
        kofzip(z) = exp(-im*p/(1-p)*π)*(z-p)^(p/(1-p))*(z-(p-1))
        kof_prime(z) = exp(-im*p/(1-p)*π)*(p/(1-p)*(z-p)^(p/(1-p)-1)*(z-(p-1)) + (z-p)^(p/(1-p)))
        u = exp(-im*p/(1-p)*π)*v^(1/(1-p))
        return newton(z->kofzip(z)-u,kof_prime,p-1+u;args...)        
    end
end

function h(z::Union{Real,Complex},c::Complex,a::Complex;args...)
    if imag(z) < 0.0
        error("h(z) is defined only for z in the upper half plane") 
    end
    p,q = reim(c)
    r,s = reim(a)
    b = (p^2*s + q^2*s - q*(r^2 + s^2))/(-q*r + p*s)
    d = a/(1-a/b)
    if isnan(z)
        return fzipinv(d,-b;args...)
    end
    return fzipinv(d,z/(1-z/b);args...)
end

function hinv(z::Union{Real,Complex},c::Complex,a::Complex;args...)
    p,q = reim(c)
    r,s = reim(a)
    b = (p^2*s + q^2*s - q*(r^2 + s^2))/(-q*r + p*s)
    d = a/(1-a/b)
    return (z->b*z/(b+z))(fzip(d,z))
end

function hcompose{T<:AbstractFloat}(A::Array{Complex{T},1},
                                    z::Union{Real,Complex};
                                    realline::Bool=false)
    zpr = z
    for i=1:2:length(A)
        realline ? zpr = real(h(zpr,A[i],A[i+1])) : zpr = h(zpr,A[i],A[i+1])
    end
    return zpr
end

function hinvcompose{T<:AbstractFloat}(A::Array{Complex{T},1},
                                       w::Union{T,Complex{T}})
    for i=length(A)-1:-2:1
        w = hinv(w,A[i],A[i+1])
    end
    return w
end

φ1(z::Union{Real,Complex},
   z0::Union{Real,Complex},
   z1::Union{Real,Complex},
   z2::Union{Real,Complex}) = im*sqrt(-((z-z2)*(z1-z0))/((z-z0)*(z1-z2)))

function φ1inv(w::Union{Real,Complex},
        z0::Union{Real,Complex},
        z1::Union{Real,Complex},
        z2::Union{Real,Complex}) 
    wsq = w*w
    return (wsq*z0*z1 + z0*z2 - wsq*z0*z2 - z1*z2)/(z0 - z1 + wsq*z1 - wsq*z2)
end

function α(z,x) 
    return π/2 + angle(abs(x/2) + im*(imag(z)/2+real(z)*(real(z)-x)/(2*imag(z))))
end

function φlast(z::Union{Real,Complex},lastpt::Real,α::Real)
    return (exp(-im*(π-α))*z/(1-z/lastpt))^(π/α) 
end

function φlastinv(w::Union{Real,Complex},lastpt::Union{Real,Complex},α::Real)
    v = w^(α/π)
    return v / (exp(-im*(π-α)) + v/lastpt)
end

#-----------------------------------------------------------------------------
# FUNCTIONS FOR DISPLAY SUPPORT  
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
        return vcat(γ,γ[1,:])
    end
end

function closepath{T<:Real}(γ::Array{T,2})
    return closepath(float(γ))
end

function plotgrid(zvals;color1=Graphics2D.blue,color2=Graphics2D.red,args...)
    lines = Graphics2D.Line[]
    for i=1:size(zvals)[1]
        push!(lines,Graphics2D.Line(hcat(real(zvals[i,:])',imag(zvals[i,:])');color=color1,args...))
    end
    for j=1:size(zvals)[2]
        push!(lines,Graphics2D.Line(hcat(real(zvals[:,j]),imag(zvals[:,j]));color=color2,args...))
    end
    return lines
end

function hyperbolictiling(f::Function;
                          rings::Integer=9,
                          rays::Integer=16,
                          rotation::Real=0.0,
                          innerradius::Real=1.0/3.0,
                          ringcolor::Array{Float64,1}=Graphics2D.blue,
                          raycolor::Array{Float64,1}=Graphics2D.red)
    points = Array{Complex64,1}[[f((1-(1-innerradius)/2^(k-1))*cos(θ+rotation) + 
                                   im*(1-(1-innerradius)/2^(k-1))*sin(θ+rotation)) 
                                 for θ=linspace(0,2π,1+rays*2^(k-1))] for k=1:rings]
    return [vcat([[Graphics2D.Line([points[i][k],points[i][k+1]];
                color=ringcolor,linesize=1-i/(rings+4)) 
            for k=1:length(points[i])-1] 
            for i=1:length(points)]...);
            vcat([[Graphics2D.Line([points[i][k],points[i+1][2*k-1]];
                color=raycolor,linesize=1-i/(rings+4)) 
            for k=1:length(points[i])-1] 
            for i=1:length(points)-1]...)]
end

function distance{T<:AbstractFloat}(p::Tuple{T,T},q::Array{T,2},r::Array{T,2})
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

function distance{T<:AbstractFloat}(p::Tuple{T,T},gamma::Array{T,2})
    m = hypot(p[1]-gamma[1,1],p[2]-gamma[1,2])
    for i=1:size(gamma)[1]-1
        d = distance(p,gamma[i,:],gamma[i+1,:])
        if d < m
            m = d
        end
    end
    return m
end

function inside{T<:AbstractFloat}(p::Tuple{T,T},gamma::Array{T,2},tol::Float64=1e-8)
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

inside{T<:AbstractFloat}(p::Tuple{T,T},gamma::Array{Complex{T},1},tol::Float64=1e-8) = 
    inside(p,hcat(reim(gamma)...),tol)

function intersect{T<:AbstractFloat}(p::Tuple{T,T},q::Tuple{T,T},r::Tuple{T,T},s::Tuple{T,T})
    a = -((p[2]*(r[1] - s[1]) + r[2]*s[1] - r[1]*s[2] + 
          p[1]*(-r[2] + s[2]))/((-p[2] + q[2])*(r[1] - s[1]) + (p[1] - 
                 q[1])*(r[2] - s[2])))
    b = (-q[2]*r[1] + p[2]*(-q[1] + r[1]) + p[1]*(q[2] - r[2]) + 
         q[1]*r[2])/((p[2] - q[2])*(r[1] - s[1]) - (p[1] - q[1])*(r[2] - s[2]))
    return 0 <= a <= 1 && 0 <= b <= 1
end

function intersect{T<:AbstractFloat}(p::Complex{T},q::Complex{T},r::Complex{T},s::Complex{T})
    # Determines whether the segment from p to q intersects the 
    # segment from r to s
    p,q,r,s = map(reim,(p,q,r,s))
    return intersect(p,q,r,s)
end

function intersect{T<:AbstractFloat}(p::Complex{T},q::Complex{T},gamma::Array{Complex{T},1})
    for i=1:length(gamma)-1
        if intersect(p,q,gamma[i],gamma[i+1])
            return true
        end
    end
    return false
end

function intersect{T<:AbstractFloat}(p::Tuple{T,T},q::Tuple{T,T},gamma::Array{T,2})
    for i=1:size(gamma)[1]-1
        if intersect(p,q,(gamma[i,1],gamma[i,2]),(gamma[i+1,1],gamma[i+1,2]))
            return true
        end
    end
    return false
end

function domain(f::Union{ConformalMap,InverseConformalMap})
    return Graphics2D.GraphicElement[Graphics2D.Point(f.center),
                                     Graphics2D.Line(closepath(f.domain))]
end

function makegrid{T<:AbstractFloat}(boundary::Array{T,2},n::Integer)
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
                       ~intersect(totalgrid[i,j],totalgrid[v...],boundary)
                    push!(lines,((i,j),v))
                end
            end
        end
    end
    return totalgrid, pointsinside, lines
end

makegrid{T<:Complex}(boundary::Array{T,1},n::Integer) = 
    makegrid(hcat(reim(boundary)...),n)

function showgrid(domain,totalgrid,pointsinside,lines,center)
    return [[Graphics2D.Point(real(center),imag(center)),
    Graphics2D.Line(closepath(domain))];
    [Graphics2D.Line([totalgrid[line[1]...],
    totalgrid[line[2]...]],linesize=0.2) for line in lines]] 
end

function showgridimage(f,totalgrid,pointsinside,lines,center)
    return [[Graphics2D.Point(0.0,0.0),Graphics2D.Circle([0.0,0.0],1)];
     [Graphics2D.Line([f(totalgrid[line[1]...]),
                       f(totalgrid[line[2]...])],linesize=0.1) for line in lines]] 
end

function visualize(CM::ConformalMap,n::Integer=40)
    totalgrid,pointsinside,lines = makegrid(closepath(CM.domain),n)
    return (showgrid(closepath(CM.domain),totalgrid,pointsinside,lines,CM.center), 
    showgridimage(p->CM(p[1] + im*p[2]),totalgrid,pointsinside,lines,CM.center))
end

function visualize(ICM::InverseConformalMap;kwargs...)
    return (hyperbolictiling(z->z;kwargs...),hyperbolictiling(z->ICM(z);kwargs...))
end


#-----------------------------------------------------------------------------


end # module
