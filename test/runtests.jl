using AMVW
using Base.Test
using Polynomials

rs = [1.0, 2, 3]
state = AMVW.amvw(poly(rs).a)
println("Mathcing?")
println(rs)
println(complex.(state.REIGS, state.IEIGS))
