## Diagonostic code
##

## make a rotator into a full matrix

function as_full{T}(a::Rotator{T}, N::Int)
    c,s = vals(a)
    i = idx(a)
    i < N || error("too big")
    A = eye(Complex{T}, N)
    A[i:i+1, i:i+1] = [c -conj(s); s conj(c)]
    A
end

function zero_out!{T}(A::Array{T}, tol=1e-12)
    A[norm.(A) .<= tol] = zero(T)
end
function zero_out!{T}(A::Array{Complex{T}}, tol=1e-12)
    for i in eachindex(A)
        c = A[i]
        cr, ci = real(c), imag(c)
        if abs(cr) < tol
            cr = zero(T)
        end
        if abs(ci) < tol
            ci = zero(T)
        end
        A[i] = complex(cr, ci)
    end
end

## diagnostic

## create Full matrix from state object. For diagnostic purposes.
# we may or may not have a diagonal matrix to keep track or
D_matrix{T}(state::ComplexRealSingleShift{T}) = diagm(state.D)
D_matrix(state::ShiftType) = I

#function Base.full{T}(state::ComplexRealSingleShift{T}, what=:A)
function Base.full{T}(state::ShiftType{T}, what=:A)
    N = state.N
    Q = as_full(state.Q[1],N+1); for i in 2:N Q = Q * as_full(state.Q[i],N+1) end
    Ct = as_full(state.Ct[1], N+1); for i in 2:N Ct =  as_full(state.Ct[i],N+1)*Ct end
    B = as_full(state.B[1],N+1); for i in 2:N B = B * as_full(state.B[i],N+1) end
    D = D_matrix(state)
    
    #    x = -vcat(state.POLY[2:state.N], state.POLY[1], 1)
    par = iseven(state.N) ? one(T) : -one(T)
    x = -vcat(state.POLY[state.N-1:-1:1], -par * state.POLY[state.N],  par * 1)
    alpha = norm(x)
    e1 = zeros(T, state.N+1); e1[1]=one(T)
    en = zeros(T, state.N+1); en[N] = one(T)
    en1 = zeros(T, state.N+1); en1[N+1] = one(T)

    rho = transpose(en1) * Ct * e1  # scalar 
    yt = -1/rho * transpose(en1)  * Ct * B * D
    # clean
    for i in eachindex(yt)
        if norm(yt[i]) < 1e-12
            yt[i] = 0
        end
    end
    
    ## we have R = Z + x = Ct * (B * D + e1 * yt)
    Z = Ct * B * D
    zero_out!(Z)
    
    x = Ct * e1 * yt 
    
    R = Z + x
    zero_out!(R)
    what == :R && return R
    
    A = Q * R
    zero_out!(A)
    A
end


# simple graphic to show march of algorithm
function show_status{T}(state::ShiftType{T})
    qs = [norm(u.s) for u in state.Q[state.ctrs.start_index:state.ctrs.stop_index]]
    minq = length(qs) > 0 ?  minimum(qs) : 0.0

    
    x = fill(".", state.N+2)
    x[state.ctrs.zero_index+1] = "α"
    x[state.ctrs.start_index+1] = "x"
    x[state.ctrs.stop_index+2] = "Δ"
    println(join(x, ""), " ($minq)")
end

## create a rotation matrix
function rotm{T}(a::T,b, i, N)
    r = eye(T, N)
    r[i:i+1, i:i+1] = [a -conj(b); b conj(a)]
    r
end

