using AMVW
const A = AMVW
using Base.Test
using Polynomials

# transformations

# givens rotation
a,b = complex(1.0, 2.0), complex(2.0, 3.0)
c,s,r = A.givensrot(a,b)
@test norm(([c -conj(s); s conj(c)] * [a,b])[2]) <= 4eps(Float64)
a,b = complex(rand(2)...), complex(rand(2)...)
c,s,r = A.givensrot(a,b)
@test norm(([c -conj(s); s conj(c)] * [a,b])[2]) <= 4eps(Float64)



# dflip
t1 = pi/3;
alpha = complex(cos(t1), sin(t1))
d = AMVW.ComplexRotator(alpha, complex(0.0, 0.0), 1)
u = one(AMVW.ComplexRotator{Float64})
AMVW.vals!(u, complex(1.0, 2.0), complex(2.0, 3.0)); AMVW.idx!(u, 2)
M = A.as_full(u, 3) * A.as_full(d,3)
A.Dflip(u, d)
M1 = A.as_full(d, 3) * A.as_full(u,3)
u = M - M1
@test  maximum(norm.(u)) <= 4eps()



##
# fuse
r1,r2 = ones(AMVW.ComplexRotator{Float64},2)
AMVW.vals!(r1, complex(1.0, 2.0), complex(2.0, 3.0)); AMVW.idx!(r1, 1)
AMVW.vals!(r2, complex(3.0, 2.0), complex(5.0, 3.0)); AMVW.idx!(r2, 1)
M = A.as_full(r1,2) * A.as_full(r2, 2)
A.fuse(r1, r2, Val{:left})
M1 = A.as_full(r1, 2)
u = M - M1
@test maximum(norm.(u)) <= 4eps()



##
# turnover
r1,r2,r3 = ones(AMVW.ComplexRotator{Float64}, 3)
AMVW.vals!(r1, complex(1.0, 2.0), complex(2.0, 3.0)); AMVW.idx!(r1, 1)
AMVW.vals!(r2, complex(3.0, 2.0), complex(5.0, 3.0)); AMVW.idx!(r2, 2)
AMVW.vals!(r3, complex(4.0, 2.0), complex(6.0, 3.0)); AMVW.idx!(r3, 1)

M = A.as_full(r1,3) * A.as_full(r2,3) * A.as_full(r3,3)
A.turnover(r1, r2, r3, Val{:right})
## we have to work a bit harder chief as we have two diagonal matrices to deal with
## 
M1 =  A.as_full(r3,3) * A.as_full(r1,3) * A.as_full(r2,3)
u = M - M1
@test maximum(norm.(u)) <= 4eps()

# reverse
A.turnover(r3, r1, r2, Val{:left})
M2 = A.as_full(r1,3) * A.as_full(r2,3) * A.as_full(r3,3)
u = M - M2
@test maximum(norm.(u)) <= 4eps()


## poly roots
# rs = [1.0, 2, 3]
# state = AMVW.amvw(poly(rs).a)
# println("Mathcing?")
# println(rs)
# println(complex.(state.REIGS, state.IEIGS))
