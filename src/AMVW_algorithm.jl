
## Main algorithm of AMV&W
## This follows that given in the paper very closely
function AMVW_algorithm{T}(state::ShiftType{T})


    it_max = 60 * state.N
    kk = 0

    while kk <= it_max

        ## finished up!
        state.ctrs.stop_index <= 0 && return
        
        check_deflation(state)
        kk += 1

##        show_status(state)

        k = state.ctrs.stop_index

        if state.ctrs.stop_index - state.ctrs.zero_index >= 2
            
            bulge_step(state)
            state.ctrs.it_count += 1
            state.ctrs.tr -= 2
            
        elseif state.ctrs.stop_index - state.ctrs.zero_index == 1

            diagonal_block(state,  k + 1)
            eigen_values(state)

            
            state.REIGS[k], state.IEIGS[k] = state.e2
            state.REIGS[k+1], state.IEIGS[k+1] = state.e1

            if k > 1
                diagonal_block(state,  k)
                eigen_values(state)
            end
            
            diagonal_block(state, 2)
            
            if state.ctrs.stop_index == 2
                diagonal_block(state, 2)
                e1 = state.A[1,1]
                state.REIGS[1] = real(e1)
                state.IEIGS[1] = imag(e1)
            end
            state.ctrs.zero_index = 0
            state.ctrs.start_index = 1
            state.ctrs.stop_index = state.ctrs.stop_index - 2
            
        elseif state.ctrs.stop_index - state.ctrs.zero_index == 0

            diagonal_block(state, state.ctrs.stop_index + 1)
            e1, e2 = state.A[1,1], state.A[2,2] 
                        
            if state.ctrs.stop_index == 1
                state.REIGS[1], state.IEIGS[1] = real(e1), imag(e1)
                state.REIGS[2], state.IEIGS[2] = real(e2), imag(e2)
                state.ctrs.stop_index = 0
            else
                state.REIGS[k+1], state.IEIGS[k+1] = real(e2), imag(e2)
                state.ctrs.zero_index = 0
                state.ctrs.start_index = 1
                state.ctrs.stop_index = k - 1
            end
        end
    end

    warn("Not all roots were found. The first $(state.ctrs.stop_index-1) are missing.")
end

## qs is [p_{n-1}, p_{n-2}, ..., p_1, p_0] for
## monic poly x^n + p_{n-1}x^{n-1} + ... + p_1 x + p_0
## returns RealDoubleShift object
function amvw{T <: Real}(qs::Vector{T})
    state = RealDoubleShift(qs)
    init_state(state)
    AMVW_algorithm(state)
    state
end

function amvw{T <: Real}(qs::Vector{Complex{T}})
    #    state = ComplexSingleShift(qs)
    state = ComplexRealSingleShift(qs)
    init_state(state)
    AMVW_algorithm(state)
    state
end

"""
Use AMVW algorithm doubleshift alorithm to find roots
of the polynomial p_0 + p_1 x + p_2 x^2 + ... + p_n x^n encoded as
`[p_0, p_1, ..., p_n]` (the same ordering used by `Polynomials`).

Returns an object of type `RealDoubleShift`.

Example: API needs work!
```
using Polynomials
x = variable()
p = poly(x - i/10 for i in 5:10)
state = amvw(p.a)
complex.(state.REIGS, state.IEIGS)
```
"""
function poly_roots{T}(ps::AbstractVector{T})
    ## roots of poly [p0, p1, ..., pn]    
    qs, k = reverse_poly(ps)

    # k is number of 0 factors
    ## simple cases
    n = length(qs)
    if n == 0 
        rts = complex.(zeros(k), zeros(k))
    elseif n == 1
        as = vcat([-real(qs[1])], zeros(k))
        bs = vcat([-imag(qs[1])], zeros(k))
        rts = complex.(as, bs)
    elseif n == 2
        if T <: Real
            b,c = -(0.5)*qs[1], qs[2]
            e1r, e1i, e2r, e2i = qdrtc(one(T), b, c)
        else
            e1r, e1i, e2r, e2i = quadratic_equation(one(T), qs[1], qs[2])
        end
        as = vcat([e1r, e2r], zeros(k))
        bs = vcat([e1i, e2i], zeros(k))
        rts = complex.(as, bs)
    else
        state = amvw(qs)
        as = vcat(state.REIGS, zeros(k))
        bs = vcat(state.IEIGS, zeros(k))
        rts = complex.(as, bs) 
    end
    return rts
end
