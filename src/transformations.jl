##
##################################################

"""
 rotations; find values
 Real Givens
 This subroutine computes c and s such that,

 [c -s] * [a, b] = [r,0]; c^2 + s^2 = 1

 and 

 r = sqrt(|a|^2 + |b|^2).

  XXX seems faster to just return r, then not
"""
function givensrot{T <: Real}(a::T,b::T)
    iszero(b) && return (sign(a) * one(T), zero(T), abs(a))
    iszero(a) && return(zero(T), -sign(b) * one(T), abs(b))

    r = hypot(a,b)
    return(a/r,-b/r,r)
end

## givens rotation
##################################################
# Compute Givens rotation zeroing b
#
# G1 [ ar + i*ai ] = [ nrm ]
# G1 [    b      ] = [     ]
#
# all variables real (nrm complex)
# returns (copmlex(cr, ci), s) with
# u=complex(cr,ci), v=s; then [u -v; v conj(u)] * [complex(ar, ai), s] has 0
#
# XXX: Could hand write this if needed, here we just use `givens` with a flip
# to get what we want, a real s part, not s part
function givensrot{T <: Real}(a::Complex{T},b::Complex{T})
    G, r = givens(b, a, 1, 2)
    G.s, -G.c, r
end
givensrot{T <: Real}(a::Complex{T},b::T) = givensrot(a, complex(b, zero(T)))


####   Operations on [,[ terms

## The zero_index and stop_index+1 point at "D" matrices
## 
## Let a D matrix be one of [1 0; 0 1] or [-1 0; 0 1] (D^2 = I). Then we have this move
## D    --->   D  (we update the rotator)
##   [       [
##
## this is `dflip`
function dflip{T}(a::RealRotator{T}, d=one(T))
    a.s = sign(d)*a.s
end

# get d from rotator which is RR(1,0) or RR(-1, 0)
function getd{T}(a::RealRotator{T})
    c, s = vals(a)
    norm(s) <= 4eps(T) || error("a is not a diagonal rotator")
    sign(c)
end


## For complex, A "D" matrix is of type Di=[alpha 0; 0 conj(alpha)] We have to
#
# D    -->    D
#   U      U 
is_diagonal{T}(r::ComplexRotator{T}) = norm(r.s) <= eps(T)
is_identity{T}(r::ComplexRotator{T}) = r.c == one(Complex{T})


# (After absorbtion we leave an I for Q, not "D" matrix as in real case)

## This is main case
#  Q           D Q
#     D --> D 
    


"""
   D  --> D
U           V
"""
function Dflip{T}(r::ComplexRotator{T}, d::ComplexRotator{T})
    is_identity(r) && error("check first for this case")
    !is_diagonal(d) && error("d must be diagonal rotator")

    # D is fixed,
    alpha = d.c
    r.s = r.s * conj(alpha)
end



## fuse combines two rotations, a and b, into one,
## XXX If using complex rotators with real sine part, need a diagonal matrix
## coming out of this

# Had switch forVal{:left} updates a, Val{:right} updates b; but this is faster
## 
function fuse{T}(a::Rotator{T}, b::Rotator{T},::Type{Val{:left}})
    #    idx(a) == idx(b) || error("can't fuse")
    u = a.c * b.c - conj(a.s) * b.s
    a.s = conj(a.c) * b.s + a.s * b.c
    a.c = u
end
function fuse{T}(a::Rotator{T}, b::Rotator{T}, ::Type{Val{:right}})
#    idx(a) == idx(b) || error("can't fuse")
    u = a.c * b.c - conj(a.s) * b.s
    b.s = conj(a.c) * b.s + a.s * b.c
    b.c = u
end


# Turnover: Q1    Q3   | x x x |      Q1
#              Q2    = | x x x | = Q3    Q2  <-- misfit=3 Q1, Q2 shift; 
#                      | x x x |
#
# misfit is Val{:right} for <-- (right to left turnover), Val{:left} for -->
#
# This is the key computation once matrices are written as rotators
## XXX we wrote this for complex rotators where sine part may be complex
# can trim down a bit if not the case, as we don't need to fuss with alpha, beta below
##function _turnover{T}(Q1::ComplexRotator{T}, Q2::ComplexRotator{T}, Q3::ComplexRotator{T})
function _turnover{T}(Q1::Rotator{T}, Q2::Rotator{T}, Q3::Rotator{T})    
#    i,j,k = idx(Q1), idx(Q2), idx(Q3)
#    (i == k) || error("Need to have a turnover up down up or down up down: have i=$#i, j=$j, k=$k")
#    abs(j-i) == 1 || error("Need to have |i-j| == 1")
    
    c1,s1 = vals(Q1)
    c2,s2 = vals(Q2)
    c3,s3 = vals(Q3)

    # key is to find U1,U2,U3 with U2'*U1'*U3' * (Q1*Q2*Q3) = I
    # do so by three Givens rotations to make (Q1*Q2*Q3) upper triangular
    
    # initialize c4 and s4
    a = conj(c1)*c2*s3 + s1*c3 
    b = s2*s3
    # check norm([a,b]) \approx 1    
    c4, s4, temp = givensrot(a,b)#, Val{true})

    # initialize c5 and s5

    a = c1*c3 - conj(s1)*c2*s3
    b = temp
    # check norm([a,b]) \approx 1
    c5, s5, alpha = givensrot(a, b)

    alpha = alpha/norm(alpha)
    c5 *= conj(alpha) # make diagonal elements 1
    s5 *= alpha
    
    # second column
    u = -c1*conj(s3) - conj(s1)*c2*conj(c3)
    v = conj(c1)*c2*conj(c3) - s1*conj(s3)
    w = s2 * conj(c3)

    a = c4*conj(c5)*v - conj(s4)*conj(c5)*w + s5*u
    b = conj(c4)*w + s4*v

    c6, s6, beta = givensrot(a,b)

    beta = beta/norm(beta)
    c6 *= conj(beta) # make diagonal elements 1
    s6 *= beta
    
    (c4, s4, c5, s5, c6, s6)
end




function turnover{T}(Q1::Rotator{T}, Q2::Rotator{T}, Q3::Rotator{T},
                     ::Type{Val{:right}})

    c4,s4,c5,s5,c6,s6 = _turnover(Q1,Q2,Q3)
    vals!(Q3, conj(c4), -s4)
    vals!(Q1, conj(c5), -s5)
    vals!(Q2, conj(c6), -s6)
    idx!(Q3, idx(Q2))  # misfit is right one
end

turnover{T}(Q1::Rotator{T}, Q2::Rotator{T}, Q3::Rotator{T}) = turnover(Q1, Q2, Q3, Val{:right})

function turnover{T}(Q1::Rotator{T}, Q2::Rotator{T}, Q3::Rotator{T},
                     ::Type{Val{:left}})
    
    c4,s4,c5,s5,c6,s6 = _turnover(Q1,Q2,Q3)
    
    vals!(Q2, conj(c4), -s4)
    vals!(Q3, conj(c5), -s5)
    vals!(Q1, conj(c6), -s6)
    idx!(Q1, idx(Q2))   # misfit is left one
end



