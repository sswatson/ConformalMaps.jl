module ConformalMaps
using Printf
# This package provides a ConfomalMap type which numerically approximates the
# conformal map from a Jordan domain specified by a list of points on the
# boundary to a disk or half-plane.

# The algorithm used is based on the paper Convergence of the Zipper algorithm
# for conformal mapping by Don Marshall and Steffen Rohde. The function names
# in this module are modeled after the notation in the paper.
# http://arxiv.org/abs/math/0605532

#-----------------------------------------------------------------------------

import Requires,
       Base.show,
       Base.inv

export ConformalMap,
       InverseConformalMap,
       domain

#-----------------------------------------------------------------------------
# Conformal Map Types
#-----------------------------------------------------------------------------

struct ConformalMap{T<:Real} # A map from the domain to the disk
    domain::Array{Complex{T},1} # Points tracing out the boundary
    data::Array{Complex{T},1} # An array (ζ in the paper) encoding the conformal map
    center::Complex{T} # The point in the domain which maps to the origin
end

struct InverseConformalMap{T<:Real} # A map from the disk to the domain
    domain::Array{Complex{T},1}
    data::Array{Complex{T},1}
    center::Complex{T}
end

function ConformalMap(domain::Array{<:Complex,1},
                      center::Complex;
                      pointspacing::Real=diameter(domain)/250)

    if domain[1] ≈ domain[end] 
        pop!(domain) 
    end
    densepoints = densify(domain,pointspacing)

    T = float(real(eltype(densepoints)))

    ζ = Complex{T}[]

    for k=3:length(densepoints)
        push!(ζ,fcompose(ζ,phi2(densepoints[k],
                                densepoints[1],densepoints[2])))
    end

    push!(ζ,fcompose(ζ,phi2(densepoints[1],densepoints[1],densepoints[2])))

    push!(ζ,densepoints[1])
    push!(ζ,densepoints[2])

    push!(ζ,philast(fcompose(ζ[1:end-3],phi2(center,ζ[end-1],ζ[end])),ζ[end-2]))

    return ConformalMap(densepoints,ζ,convert(Complex{T},center))
end

function ConformalMap(domain::Array{<:Complex,1},
                      center::Real;
                      kwargs...)
    return ConformalMap(domain,complex(center);kwargs...)
end

function ConformalMap(domain::Array{<:Real,2},
                      center::Union{Real,Complex};
                      kwargs...)
    return ConformalMap(domain[:,1] + im*domain[:,2],complex(center);kwargs...)
end

function ConformalMap(domain::Array{Tuple{<:Real,<:Real},1},
                      center::Union{Real,Complex};
                      kwargs...)
    complexdomain = [a for (a,b) in domain] + im*[b for (a,b) in domain]
    return ConformalMap(complexdomain,center;kwargs...)
end

function InverseConformalMap(domain::Array{Complex{T},1},
                       center::Complex{T};kwargs...) where T<:Real
    return inv(ConformalMap(domain,center;kwargs...))
end

function inv(CM::ConformalMap)
    return InverseConformalMap(CM.domain,CM.data,CM.center)
end

function inv(ICM::InverseConformalMap)
    return ConformalMap(ICM.domain,ICM.data,ICM.center)
end

function (CM::ConformalMap)(z::Union{Real,Complex})
    ζ = CM.data
    return H_to_D((philast(fcompose(ζ[1:end-4],
                                    phi2(z,ζ[end-2],ζ[end-1])),ζ[end-3])),ζ[end])
end

function (ICM::InverseConformalMap)(z::Union{Real,Complex})
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
        return string(shortformat(real(z)),
                      signbit(imag(z)) ? "-" : "+",
                      shortformat(abs(imag(z))),"im")
    end
end

function show(io::IO,CM::Union{ConformalMap,InverseConformalMap})
    if typeof(CM) <: ConformalMap
        print(io,"ConformalMap{")
    else
        print(io,"InverseConformalMap{")
    end
    print(io,typeof(CM).parameters[1])
    print(io,"}(\n")
    print(io," "^3,"[",shortformat(CM.domain[1]))
    print(io,",\n")
    print(io," "^4*shortformat(CM.domain[2]))
    print(io,",\n"," "^4,"⋮\n")
    print(io," "^4*shortformat(CM.domain[end-1]))
    print(io,",\n")
    print(io," "^4*shortformat(CM.domain[end]))
    print(io,"],\n"," "^4,"center=")
    print(io,shortformat(CM.center))
    print(io,")")
end

#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# COMPLEX NUMBER CONVERSION
#-----------------------------------------------------------------------------
# for conveniently splitting complex numbers and arrays of them into
# real and imaginary parts
realify(A::Array{T,1}) where {T<:Complex} = hcat(real(A),imag(A))
realify(z::Complex{T}) where {T<:AbstractFloat} = (real(z),imag(z))
realify(x::Real) = (x,0.0)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# EXPLICIT CONFORMAL MAPS
#-----------------------------------------------------------------------------
# f and finv are the basic slit-domain conformal maps composed to
# construct the global conformal map
function f(z::Union{T,Complex{T}},
           a::Union{T,Complex{T}};
           tol::Float64=1e-12) where T<:AbstractFloat
    if isnan(z)
        return -im*sign(imag(a))*abs(a)^2*sqrt(1/real(a)^2+1/imag(a)^2)
    elseif abs(imag(a)) < tol
        return sqrt(z^2 - a^2)
    else
        return sqrt((z/(1-z/(im*abs(a)^2/imag(a))))^2-abs(a)^4/real(a)^2)
    end
end

function finv(w::Union{T,Complex{T}},
              a::Union{T,Complex{T}};
              tol::Float64=1e-12) where T<:AbstractFloat
    if abs(imag(a)) < tol
        return sqrt(w^2+a^2)
    elseif isnan(w)
        return im*abs(a)^2/imag(a)
    else
        return abs(a)^2/imag(a)*sqrt(w^2+(abs(a)^2/real(a))^2) /
               (abs(a)^2/imag(a)-im*sqrt(w^2+(abs(a)^2/real(a))^2))
    end
end

function fcompose(A::Array{<:Complex,1},
                  z::Union{Real,Complex})
    for i = 1:length(A)
        z = f(z,A[i])
    end
    return z
end

function finvcompose(A::Array{<:Complex,1},
                     w::Union{Real,Complex})
    for i = length(A):-1:1
        w = finv(w,A[i])
    end
    return w
end

# The first and last maps (phi2 and philast) in Figure 3 in [MR06] have special forms

function phi2(z::Union{Real,Complex},
              point1::Union{Real,Complex},
              point2::Union{Real,Complex})
    if z == point1
        return convert(typeof(z),NaN)
    else
        return sqrt((z-point2)/(z-point1))
    end
end

phi2inv(w::Union{Real,Complex},
        point1::Union{Real,Complex},
        point2::Union{Real,Complex}) = (w^2*point1-point2)/(w^2-1.0)
philast(z::Union{Real,Complex},
        ζlast::Union{Real,Complex}) = -(z/(1-z/(ζlast)))^2
philastinv(w::Union{Real,Complex},
           ζlast::Union{Real,Complex}) = (w*ζlast+im*sqrt(w)*ζlast^2)/(w+ζlast^2)

# Explicit conformal maps between the upper half H and the disk D
H_to_D(z::Union{Real,Complex},
       a::Union{Real,Complex}) = (z-a)/(z-conj(a))

D_to_H(w::Union{Real,Complex},
       a::Union{Real,Complex}) = (w*conj(a)-a)/(w-1);
#-----------------------------------------------------------------------------



#-----------------------------------------------------------------------------
# GEODESIC ALGORITHM CONFORMAL MAPS
#-----------------------------------------------------------------------------

"""
Insert additional points along each segment of the polygonal domain boundary
"""
function densify(A::Array{<:Complex,1},pointspacing::Real)
    densearray = float(eltype(A))[]
    L = length(A)
    for j = 1:length(A)
        n = ceil(Integer,abs(A[mod1(j,L)]-A[mod1(j+1,L)])/pointspacing)
        for k=1:n
            push!(densearray,(n+1-k)/n * A[mod1(j,L)] + (k-1)/n * A[mod1(j+1,L)])
        end
    end
    return densearray
end

"""
Return the largest distance between any pair of points in `points` 
"""
function diameter(points)
    m = 0.0
    for p in points
        for q in points
            m = max(m,hypot(p,q))         
        end
    end
    return m
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
        kofzip2(z) = exp(-im*p/(1-p)*π)*(z-p)^(p/(1-p))*(z-(p-1))
        kof_prime2(z) = exp(-im*p/(1-p)*π)*(p/(1-p)*(z-p)^(p/(1-p)-1)*(z-(p-1)) + (z-p)^(p/(1-p)))
        u = exp(-im*p/(1-p)*π)*v^(1/(1-p))
        return newton(z->kofzip2(z)-u,kof_prime2,p-1+u;args...)
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

function hcompose(A::Array{Complex{T},1},
                  z::Union{Real,Complex};
                  realline::Bool=false) where T<:AbstractFloat
    zpr = z
    for i=1:2:length(A)
        realline ? zpr = real(h(zpr,A[i],A[i+1])) : zpr = h(zpr,A[i],A[i+1])
    end
    return zpr
end

function hinvcompose(A::Array{Complex{T},1},
                     w::Union{T,Complex{T}}) where T<:AbstractFloat
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

Requires.@require AsyPlots="77e5a97a-5ef9-58df-9d21-21957d92d960" begin
    export visualize
    export combine
    include("visualization.jl")
end

end # module
