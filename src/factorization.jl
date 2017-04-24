# ## Factor R into Ct, B, Dn terms
# # R = 1 0 0 -v2 0
# #     0 1 0 -v3 0
# #     0 0 1 -v1 1
# #     0 0 0   0 0
# # So we represent R by [v1, v2, v3]
# # factors R so that D * Ct * B = Z
# function R1_factorization{T}(ps::Vector{Complex{T}}, Ct, B, beta=one(T))
#     N = length(ps)
#     par = iseven(N) ? one(Complex{T}) : -one(Complex{T})
    
#     ## Working, but not quite what is in DFCC code
#     ## there -par*ps[N], par*one(T); C is -conj(c), -s
#     ## B[N] = par, -par...
#     ## Here we use: C' * B * D = Z (and not D C' B = Z)
#     c, s, temp = givensrot(par * conj(ps[N]), -par * one(Complex{T}))

#     nrm = norm(c)
#     alpha = c/nrm

#     println("Diagonal: alpha=$alpha, beta=$beta")
    
#     alpha = alpha * beta
    
#     vals!(Ct[N], conj(c), -s);
#     idx!(Ct[N], N)

#     #    vals!(B[N], -par*s, par*conj(c))
#     # B * Dn * Dn'
#     vals!(B[N], -par*s*alpha, par*norm(c))
#     idx!(B[N], N)
    
#     for ii in 2:N
#         c, s, temp = givensrot(-ps[ii-1], temp)
#         vals!(Ct[N-ii + 1], conj(c*alpha), -s)
#         idx!(Ct[N-ii+1], N-ii+1)
        
#         vals!(B[N-ii + 1], c*alpha, s)
#         idx!(B[N-ii+1], N-ii+1)
#     end

#     alpha
# end

# ## No D passed, we do a Ct * B factorization without a rotation to keep Complex Real
# ## this can be Real/Real or Complex/Complex
# function R1_factorization{T <: Real}(ps::Vector{T}, Ct, B)
#     N = length(ps)
#     par = iseven(N) ? one(T) : -one(T)
    
#     ## Working, but not quite what is in DFCC code
#     ## there -par*ps[N], par*one(T); C is -conj(c), -s
#     ## B[N] = par, -par...
#     ## Here we use: C' * B * D = Z (and not D C' B = Z)
#     c, s, temp = givensrot(par * ps[N], -par * one(T))

    
#     vals!(Ct[N], c, -s);
#     idx!(Ct[N], N)

#     #    vals!(B[N], -par*s, par*conj(c))
#     # B * Dn * Dn'
#     vals!(B[N], -par * s, par * c)
#     idx!(B[N], N)
    
#     for ii in 2:N
#         c, s, temp = givensrot(-ps[ii-1], temp)
#         vals!(Ct[N-ii + 1], c, -s)
#         idx!(Ct[N-ii+1], N-ii+1)
        
#         vals!(B[N-ii + 1], c, s)
#         idx!(B[N-ii+1], N-ii+1)
#     end

#     one(T) # alpha
# end

## Factor
##
##  1 0 -p1                   1 0 0 0 0       -p1
##  0 1 -p2  --> Yn + X -->   0 1 0 0 0       -p2
##  0 0 -p3                   0 0 1 0 -1  +   -p3
##                            0 0 0 1 0       - 1
## we pass in ps = [-p1, -p2, ..., -pn, -1] the whol column including 1
## We can leave as Ct * B * D or D * Ct * B (the defaul)
function R_factorization{T}(ps::Vector{Complex{T}}, Ct, B, side=Val{:left})

    N = length(ps) - 1 # ps has 1 in it?
    c,s,tmp = givensrot(conj(ps[N]), -one(Complex{T}))

    nrm = norm(c)
    alpha = c/nrm
    
    alpha = alpha
    
    vals!(Ct[N], conj(c), -s);
    idx!(Ct[N], N)
    
    vals!(B[N], -s*alpha, norm(c))
    idx!(B[N], N)

    ## do we leave on left? If so, it must pass through Cts and B
    gamma = side == Val{:left} ? alpha : one(Complex{T})

    for i in (N-1):-1:1
        c,s,tmp = givensrot(ps[i], tmp)
        vals!(Ct[i], conj(c*gamma), -s)
        idx!(Ct[i], i)

        vals!(B[i], c*gamma, s)
        idx!(B[i], i)
    end

    side == Val{:left} ? alpha : conj(alpha)
end


## Factor R =
# 1 0 0 p1     1 0 0 0       p1
# 0 1 0 p2 --> 0 1 0 0  -->  p2
# 0 0 1 p3     0 0 0 -1      p3
#              0 0 1 0       -1
function R_factorization{T <: Real}(ps::Vector{T}, Ct, B, side=Val{:not_applicable})
    N = length(ps) - 1
    c,s,tmp = givensrot(ps[N], - one(T))

    vals!(Ct[N], c, -s);
    idx!(Ct[N], N)
    vals!(B[N], -s, c)
    idx!(B[N], N)
    
   
    for i in (N-1):-1:1
        c,s,tmp = givensrot(ps[i], tmp)
        vals!(Ct[i], c, -s)
        idx!(Ct[i], i)

        vals!(B[i], c, s)
        idx!(B[i], i)
    end
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
    #XXX    alpha =  R_factorization(ps, state.Ct, state.B)
    
    alpha =  R_factorization(ps, state.Ct, state.B)    
    alpha
end

# Here we have V = D Ct B; W = D1 Ct1 B1
function init_triu{T, St, Tw}(state::FactorizationType{T, St, Val{:HasPencil}, Tw}, decompose)

    ps, qs = decompose(state.POLY)
    beta =  R_factorization(qs, state.Ct1, state.B1)
    
    ## if SingleShift, we need to stash beta into D1:
    if St == Val{:SingleShift}
        state.D1[state.N] = beta
        state.D1[state.N+1] = conj(beta)
    end
    
    alpha =  R_factorization(ps, state.Ct, state.B)
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




