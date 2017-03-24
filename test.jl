
## Drivers...

using Polynomials
## ps a poly with coefficients p0 + p1x +p2x^2 + p3x^3 + p4x^4 [p0, p1, ..., p4]
damvw(p::Poly) = damvw(p.a)

function poly_roots(p::Poly)
    state = damvw(p.a)
    AMVW(state)
    poly_roots(state)
end


using Base.Profile
using Polynomials
T = Float64
#T = BigFloat
x = variable(T)
p = prod(x - i/10 for i in 1:20)
#p =  prod(x^2 + i for i in 1:5)
#p = poly(linspace(.1,1,20))
state = damvw(p)
AMVW.AMVW(state)



# ## warmed up

n = 5
p = poly(linspace(.5,1,n))
println(n)
state = damvw(p)
@time DAMVW.AMVW(state)

Profile.clear()
p = poly(linspace(.5,1,n))
state = damvw(p)
@profile DAMVW.AMVW(state)
Profile.print(format=:flat, sortedby=:count)

# n = 10
# as = zeros(5n)
# for i in 1:5n
#     println("doing $i")
#     p = poly(linspace(.2, 1.0, i+2))
#     state = damvw(p)
#     a = time(); DAMVW.AMVW(state); b = time() - a
#     as[i] = b
# end
