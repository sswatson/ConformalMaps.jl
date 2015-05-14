module ConformalMaps

# This package provides an object which numerically approximates the conformal map 
# from a Jordan domain specified by a list of points on the boundary to a disk or 
# half-plane. 

# The algorithm used is based on the paper Convergence of the Zipper algorithm 
# for conformal mapping by Don Marshall and Steffen Rohde. 
# http://arxiv.org/abs/math/0605532

#-----------------------------------------------------------------------------
# BASIC USAGE 
#-----------------------------------------------------------------------------

## Computing the list of points ζ can take a long time for a complicated domain
# ζ = initialize_conformal_map(Float64[1 0; 0 1; -1 0; 0 -1];resolution=100)

## Once ζ computed, function calls are fast. 
# f(z) = conformalmap(ζ,z,0.0)

## To display the domain with curves depicting the conformal map  
# domain = Complex{Float64}[invconformalmap(ζ,r*exp(im*θ),center) 
#                  for r=linspace(0,0.99999,100),θ=linspace(0,2π,100)] 
# showgraphics([plotgrid(domain),Line(closepath(mypoints))])
#-----------------------------------------------------------------------------


import Graphics2D

export initializeconformalmap,
       conformalmap,
       invconformalmap,
       plotgrid, 
       closepath,
       hyperbolictiling


#-----------------------------------------------------------------------------
# COMPLEX NUMBER CONVERSION
#-----------------------------------------------------------------------------
# for conveniently splitting complex numbers and arrays of them into 
# real and imaginary parts 
realify{T<:Complex}(A::Array{T,1}) = hcat(real(A),imag(A))
realify{T<:FloatingPoint}(z::Complex{T}) = (real(z),imag(z))
realify(x::Real) = (x,0.0)
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
# EXPLICIT CONFORMAL MAPS
#-----------------------------------------------------------------------------
# f and finv are the basic slit-domain conformal maps composed to 
# construct the global conformal map 
function f{T<:FloatingPoint}(z::Union(T,Complex{T}),
                             a::Union(T,Complex{T});
                             tol::Float64=1e-12)
    if isnan(z)
        return -im*sign(imag(a))*abs(a)^2*sqrt(1/real(a)^2+1/imag(a)^2)
    elseif abs(imag(a)) < tol
        return sqrt(z^2 - a^2)
    else
        return sqrt((z/(1-z/(im*abs(a)^2/imag(a))))^2-abs(a)^4/real(a)^2)
    end
end

function finv{T<:FloatingPoint}(w::Union(T,Complex{T}),
                                a::Union(T,Complex{T});
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

function fcompose{T<:FloatingPoint}(A::Array{Complex{T},1},
                                    z::Union(FloatingPoint,Complex))
    for i=1:length(A)
        z = f(z,A[i])
    end
    return z
end

function finvcompose{T<:FloatingPoint}(A::Array{Complex{T},1},
                                       w::Union(T,Complex{T}))
    for i=length(A):-1:1
        w = finv(w,A[i])
    end
    return w
end

# The first and last maps (phi2 and philast) in Figure 3 in [MR06] have special forms

function phi2(z::Union(Complex,FloatingPoint),
              point1::Union(Complex,FloatingPoint),
              point2::Union(Complex,FloatingPoint))
    if z == point1
        return convert(typeof(z),NaN)
    else
        return sqrt((z-point2)/(z-point1))
    end
end

phi2inv(w::Union(Complex,FloatingPoint),
        point1::Union(Complex,FloatingPoint),
        point2::Union(Complex,FloatingPoint)) = (w^2*point1-point2)/(w^2-1.0)
philast(z::Union(Complex,FloatingPoint),
        ζlast::Union(Complex,FloatingPoint)) = -(z/(1-z/(ζlast)))^2
philastinv(w::Union(Complex,FloatingPoint),
           ζlast::Union(Complex,FloatingPoint)) = (w*ζlast+im*sqrt(w)*ζlast^2)/(w+ζlast^2)

# Explicit conformal maps between the upper half H and the disk D
H_to_D(z::Union(Complex,FloatingPoint),
       a::Union(Complex,FloatingPoint)) = (z-a)/(z-conj(a))

D_to_H(w::Union(Complex,FloatingPoint),
       a::Union(Complex,FloatingPoint)) = (w*conj(a)-a)/(w-1);
#-----------------------------------------------------------------------------




#-----------------------------------------------------------------------------
# ZIPPER ALGORITHM CONFORMAL MAPS
#-----------------------------------------------------------------------------

# This function adds additional points along the segments defining 
# the boundary of the domain
function densify{T<:FloatingPoint}(A::Array{T,2},n::Integer)
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

function initializeconformalmap{T<:FloatingPoint}(points::Array{T,2};resolution=1)
    densepoints = densify(points,resolution)
    complexpoints = densepoints[:,1] + im*densepoints[:,2];
    
    ζ = Complex{typeof(points[1,1])}[]

    for k=3:length(complexpoints)
        push!(ζ,fcompose(ζ,phi2(complexpoints[k],complexpoints[1],complexpoints[2])))
    end

    push!(ζ,fcompose(ζ,phi2(complexpoints[1],complexpoints[1],complexpoints[2])))

    push!(ζ,complexpoints[1])
    push!(ζ,complexpoints[2])
    
    return ζ
end

function initializeconformalmap{T<:Real}(points::Array{T,2};args...)
    return initializeconformalmap(float(points);args...)
end

function initializeconformalmap{T<:Real}(points::Array{Complex{T},1};args...) 
    return initializeconformalmap(hcat(real(points),imag(points));args...)
end

function conformalmap{T<:FloatingPoint}(ζ::Array{Complex{T},1},
                                        z::Union(T,Complex{T}),
                                        center::Union(T,Complex{T}))
    a = philast(fcompose(ζ[1:end-3],phi2(center,ζ[end-1],ζ[end])),ζ[end-2])
    return (x->(x-a)/(x-conj(a)))(philast(fcompose(ζ[1:end-1],phi2(z,ζ[end-1],ζ[end])),ζ[end]))
end

conformalmap{T<:FloatingPoint}(ζ::Array{Complex{T},1},
                               p::(T,T),
                               a::Union(T,Complex{T})) = conformalmap(ζ,p[1] + im*p[2],a)

function invconformalmap{T<:FloatingPoint}(ζ::Array{Complex{T},1},
                         w::Union(Complex,FloatingPoint),
                         center::Union(Complex,FloatingPoint))
    a = philast(fcompose(ζ[1:end-3],phi2(center,ζ[end-1],ζ[end])),ζ[end-2])
    return phi2inv(finvcompose(ζ[1:end-3],philastinv(D_to_H(w,a),ζ[end-2])),ζ[end-1],ζ[end])
end
#-----------------------------------------------------------------------------


#-----------------------------------------------------------------------------
# FUNCTIONS FOR DISPLAY SUPPORT  
#-----------------------------------------------------------------------------

function closepath{T<:FloatingPoint}(γ::Array{Complex{T},1})
    return closepath(hcat([real(a) for a in γ], [imag(a) for a in γ]))
end

function closepath{T<:Real}(γ::Array{Complex{T},1})
    return closepath(float(γ))
end

function closepath{T<:FloatingPoint}(γ::Array{T,2})
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
                           rings::Integer=8,
                           rays::Integer=2,
                           rotation::Real=0.0,
                           innerradius::Real=0.5,
                           ringcolor::Array{Float64,1}=Graphics2D.blue,
                           raycolor::Array{Float64,1}=Graphics2D.red)
    points = Array{Complex64,1}[[f((1-(1-innerradius)/2^(k-1))*cos(θ+rotation) + 
            im*(1-(1-innerradius)/2^(k-1))*sin(θ+rotation)) 
            for θ=linspace(0,2π,1+rays*2^(k-1))] for k=1:rings]
    return [vcat([[Graphics2D.Line([points[i][k],points[i][k+1]];
                color=ringcolor,linesize=1-i/(rings+4)) 
            for k=1:length(points[i])-1] 
            for i=1:length(points)]...),
            vcat([[Graphics2D.Line([points[i][k],points[i+1][2*k-1]];
                color=raycolor,linesize=1-i/(rings+4)) 
            for k=1:length(points[i])-1] 
            for i=1:length(points)-1]...)]
end


#-----------------------------------------------------------------------------


end # module
