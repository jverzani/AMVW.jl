## Factorization code

## Factor
##
##  1 0 -p1                   1 0 0 0 0       -p1
##  0 1 -p2  --> Yn + X -->   0 1 0 0 0       -p2
##  0 0 -p3                   0 0 1 0 -1  +   -p3
##                            0 0 0 1 0       - 1
## we pass in ps = [-p1, -p2, ..., -pn, -1] the whol column including 1
## We can leave as Ct * B * D or D * Ct * B (the default)
##
## XXX -- user should not have to worry about use of `par` and `par * one(T)` for last two terms
## should do this here, not in the `decompose` functions
##
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




