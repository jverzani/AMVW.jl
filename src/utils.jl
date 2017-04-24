function deflate_leading_zeros{T}(ps::Vector{T})
    ## trim any 0s from the end of ps
    N = findlast(!iszero, ps)
    K = findfirst(!iszero, ps)

    N == 0 && return(zeros(T,0), length(ps))
    ps = ps[K:N]
    ps, K-1
end

## take poly [p0, p1, ..., pn] and return
## [q_m-1, q_m-2, ..., q0], k
## where we trim of k roots of 0, and then make p monic, then reverese
## monomial x^5
function reverse_poly{T}(ps::Vector{T})
    # assume we have called deflate_leading_zeros
    qs = reverse(ps./ps[end])[2:end]
    qs
end

#
function quadratic_equation{T <: Real}(a::T, b::T, c::T)   
    qdrtc(a, -(0.5)*b, c)
end

## make more robust
function quadratic_equation{T}(a::Complex{T}, b::Complex{T}, c::Complex{T})
    d = sqrt(b^2 - 4*a*c)
    e1 = (-b + d)/(2a); e2 = (-b-d)/(2a)
    return (real(e1), imag(e1), real(e2), imag(e2))
    
end

## Kahan quadratic equation with fma
##  https://people.eecs.berkeley.edu/~wkahan/Qdrtcs.pdf

## solve ax^2 - 2bx + c
function qdrtc{T <: Real}(a::T, b::T, c::T)
    # z1, z2 roots of ax^2 - 2bx + c
    d = discr(a,b,c)  # (b^2 - a*c), as 2 removes 4
    
    if d <= 0
        r = b/a  # real
        s = sqrt(-d)/a #imag
        return (r,s,r,-s)
    else
        r = sqrt(d) * (sign(b) + iszero(b)) + b
        return (r/a, zero(T), c/r, zero(T))
    end
end

## more work could be done here.
function discr{T}(a::T,b::T,c::T)
    pie = 3.0 # depends on 53 or 64 bit...
    d = b*b - a*c
    e = b*b + a*c

    pie*abs(d) > e && return d

    p = b*b
    dp = muladd(b,b,-p)
    q = a*c
    dq = muladd(a,c,-q)

    (p-q) + (dp - dq)
end

##
# solve degree 2 or less case
## COMPLEX VALUSE XXX
function solve_simple_cases(state)
#    println("Simple case setting eigen value")
    if N == 0
        state.FLAG = -1
        return
    elseif N == 1
        state.FLAG = 0
        N == 1 && (state.REIGS[1] = -state.POLY[1])
        return
    elseif N == 2
        # quadratic formula
        c,b,a = state.POLY[1], state.POLY[2], 1.0

        tr = -b
        disc = b^2 - 4.0*c

        if disc < 0
            state.REIGS[1] = -b/2.0
            state.IEIGS[1] = sqrt(-disc)/2.0
            state.REIGS[2] = state.REIGS[1]
            state.IEIGS[2] = -state.IEIGS[1]
        else
            u,v = tr + sqrt(disc), tr - sqrt(disc)
            if abs(u) < abs(v)
                u,v = v, u
            end
            if u == 0
                ## nothing to do
            else
                state.REIGS[1] = u/2.0
                state.REIGS[2] = c/state.REIGS[1]
            end
        end
    end
end

