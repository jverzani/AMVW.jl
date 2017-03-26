## Bulge chasing
## RealDoubleShift
##
## There are two rotators, U, V, to chase through the matrix using the following operations
##                                   [         [
## a unitary transform: basically U    [   -->   [     U; as we just hit both sides by U' A U and U' U is I
##                                       [          [
##
##
## A fuse   [ [ -> [
##
## A turnover    [   M -->   [     where M moves through a descending or ascending structure
##                 [       M   [
##
## a "D" flip:  D    --->     D
##                [        [
##


## The bulge is created by  (A-rho1) * (A - rho2) * e_1 where rho1 and rho2 are eigenvalue or random
function create_bulge{T}(state::RealDoubleShift{T})

    if mod(state.ctrs.it_count, 15) == 0
        
        t = rand() * pi
        re1, ie1 = cos(t), sin(t)
        re2, ie2 = re1, -ie1
        
        vals!(state.U, re1, ie1); idx!(state.U, state.ctrs.start_index)
        vals!(state.Ut, re1, -ie1); idx!(state.Ut, state.ctrs.start_index)
        
        vals!(state.V, re2, ie2); idx!(state.V, state.ctrs.start_index + 1)
        vals!(state.Vt, re2, -ie2); idx!(state.Vt, state.ctrs.start_index + 1)        
        
    else

        # compute (A-rho1) * (A - rho2) * e_1
        # find e1, e2

        flag = diagonal_block(state, state.ctrs.stop_index+1)
        eigen_values(state)        
        l1r, l1i = state.e1
        l2r, l2i =  state.e2

        # find first part of A[1:3, 1:2]
        Bk = state.Bk 
        flag = flag | diagonal_block(state,  state.ctrs.start_index+1)

        bk11, bk12 = state.A[1,1], state.A[1,2]
        bk21, bk22 = state.A[2,1], state.A[2,2]

        # find last part
        flag = flag | diagonal_block(state, state.ctrs.start_index+2)
        #        Bk[3,2] = state.A[2, 1]
        bk32 = state.A[2,1]

        # an issue...
        # if isnan(l1r) || isnan(l1i) || isnan(l2r) || isnan(l2i)
        #     ## eigvals gone awry
        #     restart(state)
        #     return create_bulge(state)
        # end

        
        
#        if !flag  # flag is false if there is an issue
#            restart(state)
#            return create_bulge(state)
#        end
        
        # make first three elements of c1,c2,c3
        # c1 = real(-l1i⋅l2i + ⅈ⋅l1i⋅l2r - ⅈ⋅l1i⋅t₁₁ + ⅈ⋅l1r⋅l2i + l1r⋅l2r - l1r⋅t₁₁ - ⅈ⋅l2i⋅t₁₁ - l2r⋅t₁₁ + t₁₁^2  + t₁₂⋅t₂₁)
        # c2 = real(-ⅈ⋅l1i⋅t₂₁ - l1r⋅t₂₁ - ⅈ⋅l2i⋅t₂₁ - l2r⋅t₂₁ + t₁₁⋅t₂₁ + t₂₁⋅t₂₂)
        # c3 = real(t₂₁⋅t₃₂)
        
        c1 = -l1i * l2i + l1r*l2r -l1r*bk11 -l2r * bk11 + bk11^2 + bk12 * bk21
        c2 = -l1r * bk21 - l2r * bk21 + bk11* bk21 + bk21 * bk22
        c3 = bk21 * bk32

        
        c, s, nrm = givensrot(c2, c3)
        j = state.ctrs.start_index + 1

        vals!(state.V, c, -s)
        idx!(state.V, j)

        c, s, tmp = givensrot(c1, nrm)

        vals!(state.U, c, -s)
        idx!(state.U, j-1)
    end

end

## make W on left side
#
# initial            Q0
# we do turnover U1'     Q1     -->    U1'        -->    U1'          -->    Q1
#                    V1'    Q2      Q1     V1' Q2     Q1     (V1'Q2)       W1  Q2
#
# With this, W will be kept on the left side until the last step, U,V
# move through left to right by one step, right to left by unitariness
#
#      Q0                       Q0                     Q0               Q0
#  U1'     Q1        U1'    Q1*         ->         U1            -->       U1*
#      V1'    Q3 ->     V1'         Q3      Q1**      V1'   Q3        W          (V1'Q3)
#
# Q0 is (p,0) rotator, p 1 or -1. We have
#    Q0  --> Q0
#  R             (r, pr2)
function prepare_bulge{T}(state::RealDoubleShift{T})
    
    # N = state.N
    # as_full(state.V', N+1)* as_full(state.U', N+1)* full(state) * as_full(state.U, N+1) * as_full(state.V, N+1) |> eigvals |> println


    k = state.ctrs.start_index

    vals!(state.Ut, state.U.c, -state.U.s); idx!(state.Ut, idx(state.U))
    vals!(state.Vt, state.V.c, -state.V.s); idx!(state.Vt, idx(state.V))

    copy!(state.W, state.Q[k])
    p = k == 1 ? one(T) : state.Q[k-1].c  #  zero index implies Q0 = RR(1,0) or RR(-1,0)
    dflip(state.W, p)
    
    turnover(state.Ut, state.Vt, state.W, Val{:right})
    fuse(state.Vt, state.Q[k+1], Val{:right})  # V' Q3
    dflip(state.Ut, p)
    vals!(state.Q[k], state.Ut.c, state.Ut.s) 
    
end

## Bulge chasing moves U, V fr from R to L through B then C then Q where an interaction with W allows
## a subsequent unitary operation to move U,V back to the right, one step down
## The case when Ct[i] and B[i] are identical allow a speed up.

function one_bulge_chase_shortcut{T}(state::RealDoubleShift{T})
    i = idx(state.V)
    # borrow Vt, Ut here to store a copy
    copy!(state.Vt, state.V)
    copy!(state.Ut, state.U)    
    
    turnover(state.B[i],    state.B[i+1], state.Vt, Val{:right})
    turnover(state.B[i-1],  state.B[i],   state.Ut, Val{:right})
    for k in -1:1
        a,b = vals(state.B[i+k])
        vals!(state.Ct[i+k], a, -b) # using copy!(Ct, B') is slower
    end

    turnover(state.Q[i],    state.Q[i+1], state.V, Val{:right})
    turnover(state.Q[i-1],  state.Q[i],   state.U, Val{:right})
    turnover(state.W,       state.V,      state.U, Val{:left})
    
end

function one_bulge_chase{T}(state::RealDoubleShift{T})
    i = idx(state.V)
    turnover(state.B[i],    state.B[i+1], state.V, Val{:right})
    turnover(state.Ct[i+1], state.Ct[i],  state.V, Val{:right})

    
    j = i - 1
    turnover(state.B[j],    state.B[j+1], state.U, Val{:right})
    turnover(state.Ct[j+1], state.Ct[j],  state.U, Val{:right})

    turnover(state.Q[i],    state.Q[i+1], state.V, Val{:right})
    turnover(state.Q[j],    state.Q[j+1], state.U, Val{:right})
    turnover(state.W,       state.V,      state.U, Val{:left})

end


function chase_bulge{T}(state::RealDoubleShift{T})

    # println("  begin chase at level $(state.V.i)")
    # as_full(state.W, state.N+1)* full(state) * as_full(state.V, state.N+1) * as_full(state.U, state.N+1) |> eigvals |> println
    
    # one step
    i = idx(state.V)

    ## When  i < tr  C_i = B_i. This happens in the early steps
    ## this means fewer turnovers, but at a price of more allocations
    while i < state.ctrs.stop_index # loops from start_index to stop_index - 1
        if i <= state.ctrs.tr
            one_bulge_chase_shortcut(state) 
        else
            one_bulge_chase(state)
        end

        i += 1

    end

    # println("end chase")
    # as_full(state.W, state.N+1)* full(state) * as_full(state.V, state.N+1) * as_full(state.U, state.N+1) |> eigvals |> println    
end


## Bulge is absorbed by moving V through, then U going throug two trips.
function absorb_bulge{T}(state::RealDoubleShift{T})

    # println("absorb 0")
    # as_full(state.W, state.N+1) * full(state) * as_full(state.V, state.N+1) * as_full(state.U, state.N+1) |> eigvals |> println    

    
    # first V goes through B, C then fuses with Q
    i = idx(state.V)

    turnover(state.B[i],     state.B[i+1], state.V, Val{:right})
    turnover(state.Ct[i+1],  state.Ct[i],  state.V, Val{:right})

    ## We may be fusing Q            P  --> (Q') 
    #                      RR(-1,0)               RR(-1,0)
    #
    
    p = getd(state.Q[i+1])
    dflip(state.V, p)
    fuse(state.Q[i], state.V, Val{:left}) # fuse Q*V -> Q


    # println("absorb 1")
    # as_full(state.W, state.N+1) * full(state) *  as_full(state.U, state.N+1) |> eigvals |> println        

    
    # Then bring U through B, C, and Q to fuse with W
    j = idx(state.U)
    turnover(state.B[j],     state.B[j+1], state.U)
    turnover(state.Ct[j+1],  state.Ct[j],  state.U)
    turnover(state.Q[j],     state.Q[j+1], state.U)
    fuse(state.W, state.U, Val{:right})

    # println("absorb 2")
    # as_full(state.U, state.N+1) * full(state) |> eigvals |> println    
    
    # similarity transformation, bring through then fuse with Q
    j = idx(state.U)
    turnover(state.B[j], state.B[j+1], state.U, Val{:right})
    turnover(state.Ct[j+1],  state.Ct[j], state.U)
    p = getd(state.Q[j+1])
    dflip(state.U, p)
    fuse(state.Q[j], state.U, Val{:left})

    # println("absorb final")
    # full(state) |> eigvals |> println    
end


##################################################

## ComplexSingleShift 
function create_bulge{T}(state::ComplexSingleShift{T})

    if mod(state.ctrs.it_count, 15) == 0
        
        t = rand() * pi
        if state.ray
            shift = complex(cos(t), sin(t))
        else
            shift = complex(cos(t), zero(T))
        end
        
    else
        
        flag = diagonal_block(state, state.ctrs.stop_index+1)
        if state.ray
            e1, e2 = eigen_values(state)
            shift = norm(state.A[2,2] - e1) < norm(state.A[2,2] - e2) ? e1 : e2
        else
            shift = state.A[2,2]
        end
        
    end

    flag = diagonal_block(state, state.ctrs.start_index+1)
    c,s,nrm = givensrot(state.A[1,1] - shift, state.A[2,1])

    vals!(state.U, conj(c), -s) # U is the inverse of what we just found, 
    idx!(state.U, state.ctrs.start_index)
    vals!(state.Ut, c, s)
    idx!(state.Ut, idx(state.U))
    nothing
end
        
  

##
##    D        D           D
## U'    Q -->   V Q -->     (VQ) 
##
## with D = D(alpha); U' = (u1, v1); V = (u1 conj(alpha), v1)
function prepare_bulge{T}(state::ComplexSingleShift{T})
    i = idx(state.Ut)
    if i > 1
        # if previously deflated, the prior is only diagonal
        # so may have trouble passing Ut to Q[i]
        Dflip(state.Ut, state.Q[i-1])
    end
    fuse(state.Ut, state.Q[i], Val{:right})
end       

##
function one_bulge_chase_shortcut{T}(state::ComplexSingleShift{T})
    ## XXX speed up goes here, as we only need turnover through B, not C
    ## savings are one fewer turnover, a few copies
    one_bulge_chase(state)
end
# Moving QCBU_i -> QC(BUi) -> QCUi1B -> Q(CUi1)B -> QUiCB -> Ui1 Q C B
function one_bulge_chase{T}(state::ComplexSingleShift{T})
    i = idx(state.U)
    turnover(state.B[i],    state.B[i+1], state.U, Val{:right})
    turnover(state.Ct[i+1], state.Ct[i],  state.U, Val{:right})
    turnover(state.Q[i],    state.Q[i+1], state.U, Val{:right})    
end

# Q Ct B D U -> Q Ct B (D U) -> Q Ct (B U) D -> Q (Ct U) B D ->
# (Q U) Ct B D -> U Q Ct B D  then wrap via unitary operation
function chase_bulge{T}(state::ComplexSingleShift{T})

    # one step
    i = idx(state.U)

    ## When  i < tr  C_i = B_i. This happens in the early steps
    ## this means fewer turnovers, but at a price of more allocations
    while i < state.ctrs.stop_index # loops from start_index to stop_index - 1
        if i <= state.ctrs.tr
            one_bulge_chase_shortcut(state) 
        else
            one_bulge_chase(state)
        end

        i += 1

    end

end

# We have Q Ct B D U -> Q Ct B (U D) -> Q Ct (B U) D -> Q (Ct U) B D ->
# (Q U) Ct B D -> Q Di Ct B D -> Q (Di Ct) B D -> Q Ct (Di B) D -> Q Ct B (Di D)
# Q CB(Ui) -> QUi C B -> (QUi
#  I
function absorb_bulge{T}(state::ComplexSingleShift{T})
    i = idx(state.U)

    turnover(state.B[i],    state.B[i+1], state.U, Val{:right})
    turnover(state.Ct[i+1], state.Ct[i],  state.U, Val{:right})
    i < state.N && Dflip(state.U, state.Q[i+1])
    fuse(state.Q[i], state.U, Val{:left})
end


##################################################

## chase bulge from top to bottom until final absorbtion 
function bulge_step{T}(state::ShiftType{T})
    create_bulge(state)
    prepare_bulge(state)
    chase_bulge(state)
    absorb_bulge(state)
end

