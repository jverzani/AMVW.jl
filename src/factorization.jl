## Factor R into Ct, B, Dn terms
# R = 1 0 0 -v2 0
#     0 1 0 -v3 0
#     0 0 1 -v1 1
#     0 0 0   0 0
# So we represent R by [v1, v2, v3]
# factors R so that D * Ct * B = Z
function R_factorization{T}(ps::Vector{Complex{T}}, Ct, B)
    N = length(ps)
    par = iseven(N) ? one(Complex{T}) : -one(Complex{T})
    
    ## Working, but not quite what is in DFCC code
    ## there -par*ps[N], par*one(T); C is -conj(c), -s
    ## B[N] = par, -par...
    ## Here we use: C' * B * D = Z (and not D C' B = Z)
    c, s, temp = givensrot(par * conj(ps[N]), -par * one(Complex{T}))

    nrm = norm(c)
    alpha = c/nrm
    
    vals!(Ct[N], conj(c), -s);
    idx!(Ct[N], N)

    #    vals!(B[N], -par*s, par*conj(c))
    # B * Dn * Dn'
    vals!(B[N], -par*s*alpha, par*norm(c))
    idx!(B[N], N)
    
    for ii in 2:N
        c, s, temp = givensrot(-ps[ii-1], temp)
        vals!(Ct[N-ii + 1], conj(c*alpha), -s)
        idx!(Ct[N-ii+1], N-ii+1)
        
        vals!(B[N-ii + 1], c*alpha, s)
        idx!(B[N-ii+1], N-ii+1)
    end

    alpha
end

## No D passed, we do a Ct * B factorization without a rotation to keep Complex Real
## this can be Real/Real or Complex/Complex
function R_factorization{T <: Real}(ps::Vector{T}, Ct, B)
    N = length(ps)
    par = iseven(N) ? one(T) : -one(T)
    
    ## Working, but not quite what is in DFCC code
    ## there -par*ps[N], par*one(T); C is -conj(c), -s
    ## B[N] = par, -par...
    ## Here we use: C' * B * D = Z (and not D C' B = Z)
    c, s, temp = givensrot(par * ps[N], -par * one(T))

    
    vals!(Ct[N], c, -s);
    idx!(Ct[N], N)

    #    vals!(B[N], -par*s, par*conj(c))
    # B * Dn * Dn'
    vals!(B[N], -par * s, par * c)
    idx!(B[N], N)
    
    for ii in 2:N
        c, s, temp = givensrot(-ps[ii-1], temp)
        vals!(Ct[N-ii + 1], c, -s)
        idx!(Ct[N-ii+1], N-ii+1)
        
        vals!(B[N-ii + 1], c, s)
        idx!(B[N-ii+1], N-ii+1)
    end

    one(T) # alpha
end

function Q_factorization{T, P}(state::FactorizationType{T, Val{:SingleShift}, P, Val{:NotTwisted}})
    N = state.N
    Q = state.Q
    for ii = 1:(N-1)
        vals!(Q[ii], zero(Complex{T}), one(T))
        idx!(Q[ii], ii)
    end
    vals!(Q[N], one(Complex{T}), zero(T)) #I
    idx!(Q[N], N)
end

function Q_factorization{T, P}(state::FactorizationType{T, Val{:DoubleShift}, P, Val{:NotTwisted}})
    N = state.N
    Q = state.Q
    for ii = 1:(N-1)
        vals!(Q[ii], zero(T), one(T))
        idx!(Q[ii], ii)
    end
    vals!(Q[N], one(T), zero(T))  # I
    idx!(Q[N], N)
end

# function Q_factorization{T}(state::ComplexRealSingleShift{T}, Q)
#     N = state.N
#     for ii = 1:(N-1)
#         vals!(Q[ii], zero(Complex{T}), one(T))
#         idx!(Q[ii], ii)
#     end
#     vals!(Q[N], one(Complex{T}), zero(T)) #I
#     idx!(Q[N], N)
# end

# function Q_factorization{T}(state::ShiftType{T}, Q)
#     N = state.N
#     for ii = 1:(N-1)
#         vals!(Q[ii], zero(T), one(T))
#         idx!(Q[ii], ii)
#     end
#     vals!(Q[N], one(T), zero(T))  # I
#     idx!(Q[N], N)
# end



# ## 
# ## initial factorization
# ## This is for complex real where we have a D matrix for phases
# function QDCB_factorization{T}(state::ComplexRealSingleShift{T})
#     Q_factorization(state, state.Q)
#     R_factorization(state.POLY, state.Ct, state.B, state.D)
# end

# function QCB_factorization{T}(state::ShiftType{T})
#     Q_factorization(state, state.Q)
#     R_factorization(state.POLY, state.Ct, state.B)
# end


## Init State
function init_state{T, St, P, Tw}(state::FactorizationType{T, St, P, Tw}, decompose)
    # (Q,D, sigma)
    # Ct, B [Ct1, B1]
    alpha = init_triu(state, decompose)
    init_Q(state, alpha)

end


function init_triu{T, ST, Tw}(state::FactorizationType{T, ST, Val{:NoPencil}, Tw}, decompose)
    # fill out Ct, B, pass back alpha
    ## We need to do something here with POLY!
    ps = decompose(state.POLY)
    alpha =  R_factorization(ps, state.Ct, state.B)
    alpha
end

function init_triu{T, St, Tw}(state::FactorizationType{T, St, Val{:HasPencil}, Tw}, decompose)
    # fill out Ct, B, pass back alpha
    # Need to so
    ps, qs = decompose(state.POLY)
    alpha =  R_factorization(ps, state.Ct, state.B)
    beta =  R_factorization(qs, state.Ct1i, state.B1i)
    ## we need to reverse the qs
    for i in eachindex(state.B1i)
        copy!(state.Ct1i[i], state.Ct1i[i]')
        copy!(state.B1i[i], state.B1i[i]')
    end
    ## how to combine alpha and beta? when St is SingleShift?
    alpha
end



function init_Q{T, P}(state::FactorizationType{T, Val{:SingleShift}, P, Val{:NotTwisted}}, alpha)
    # define Q, D
    Q_factorization(state)
    state.D[state.N] = alpha
    state.D[state.N+1] = conj(alpha)
end

function init_Q{T, P}(state::FactorizationType{T, Val{:DoubleShift}, P, Val{:NotTwisted}}, alpha)
    # define Q, no D
    Q_factorization(state)
end



















#init_state{T}(state::ComplexRealSingleShift{T}) = QDCB_factorization(state)
#init_state{T}(state::ShiftType{T}) = QCB_factorization(state)

## Pencils are more involved
## They have a decomposition, such as this basic one
# function simple_pencil{T}(ps::Vector{T})
#     n = length(ps) - 1
#     vs = ps[1:end-1]
#     ws = zeros(T, n)
#     ws[end] = ps[end]
#     vs, ws
# end

# ## and we initialize a bit differently
# function init_state{T}(state::ComplexRealSingleShiftPencil{T}, decompose=simple_pencil)
#     ps = state.POLY
#     n = state.N - 1
#     ## vs, ws both length n with
#     ## vs[1] = ps_0; ws[n] = ps_n; vs[i+1] + ws[i] = ps_i

#     ## Need to generalize this
#     vs, ws = decompose(ps)

#     Q_factorization(state, state.Q)
#     # do W first, as we need to move D
#     R_factorization(vs, state.Ctp, state.Bp, state.D)
#     Dn = ComplexRealRotator(state.D[n], zero(T), n)
#     R_factorization(ws, state.Ct, state.B, state.D)

#     # now we have QD(Ct*B) *(Dn*Ctp*Bp)^(-1) = QD(Ct*B) * (Ctp *Bp)^(-1) * Dn^(-1). Unitary transform to get
#     # Dn^(-1) * Q D * 
#     # Q * (Dn * D) *(Ct*B) * (Ctp * Bp)^(-1)  # move Dn past Q
#     turnover(Dn, state.Q[n-1], state.Q[n], Val{:left})
#     i = idx(Dn)
#     state.D[i] *= Dn.c
#     state.D[i+1] *= conj(Dn.c)

# end





# # If there is an issue, this function can be used to resetart the algorithm
# # could be merged with init_state?
# function restart{T}(state::ShiftType{T})
#     # try again
#     init_state(state)
    
#     for i in 1:state.N
#         state.REIGS[i] = state.IEIGS[i] = zero(T)
#     end
#     state.ctrs.zero_index = 0
#     state.ctrs.start_index = 1
#     state.ctrs.stop_index = state.N - 1
#     state.ctrs.it_count = 0
#     state.ctrs.tr = state.N - 2
# end




