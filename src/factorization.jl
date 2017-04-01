## 
## initial factorization
## This is for complex real where we have a D matrix for phases
function QDCB_factorization{T}(state::ComplexRealSingleShift{T})

    N, ps= state.N, state.POLY
    par = iseven(N) ? one(T) : -one(T)

    
    Q, Ct, B = state.Q, state.Ct, state.B
    state.D[:] = ones(Complex{T}, N+1)
    Dn = ones(Complex{T}, 2)
    
    for ii = 1:(N-1)
        vals!(Q[ii], zero(Complex{T}), one(T))
        idx!(Q[ii], ii)
    end
    vals!(Q[N], one(Complex{T}), zero(T))
    idx!(Q[N], N)


    ## Working, but not quite what is in DFCC code
    ## there -par*ps[N], par*one(T); C is -conj(c), -s
    ## B[N] = par, -par...
    c, s, temp = givensrot(par * conj(ps[N]), -par * one(T)) # <<<- conj(ps[N])!!

    nrm = norm(c)
    alpha = c/nrm
    
    vals!(Ct[N], conj(c), -s);
    idx!(Ct[N], N)

    vals!(B[N], -par*s*alpha, par*norm(c))
    idx!(B[N], N)
    
    state.D[N] = alpha
    state.D[N+1] = conj(alpha)
    
    for ii in 2:N
        c, s, temp = givensrot(-ps[ii-1], temp)
        vals!(Ct[N-ii + 1], conj(c*alpha), -s)  # note alpha
        idx!(Ct[N-ii+1], N-ii+1)
        
        vals!(B[N-ii + 1], c*alpha, s)          # note alpha
        idx!(B[N-ii+1], N-ii+1)
    end
end
init_state{T}(state::ComplexRealSingleShift{T}) = QDCB_factorization(state)

## 
## initial factorization for
## RealDoubleShift and ComplexComplex
function QCB_factorization{T}(state::ShiftType{T})

    N, ps= state.N, state.POLY
    par = iseven(N) ? one(T) : -one(T)

    
    if isa(state, RealDoubleShift)
        const ZERO, ONE = zero(T), one(T)
    else
        const ZERO, ONE = zero(Complex{T}), one(Complex{T})
    end
    

    
    Q, Ct, B = state.Q, state.Ct, state.B

    for ii = 1:(N-1)
        vals!(Q[ii], ZERO, ONE)
        idx!(Q[ii], ii)
    end
    vals!(Q[N], ONE, ZERO)
    idx!(Q[N], N)


    ## Working, but not quite what is in DFCC code
    ## there -par*ps[N], par*one(T); C is -conj(c), -s
    ## B[N] = par, -par...
    c, s, temp = givensrot(par * ps[N], -par * one(T))
    vals!(Ct[N], conj(c), -s); idx!(Ct[N], N)

    vals!(B[N], -par*s, par*conj(c))     
    idx!(B[N], N)
    
    for ii in 2:N
        c, s, temp = givensrot(-ps[ii-1], temp)
        vals!(Ct[N-ii + 1], conj(c), -s)
        idx!(Ct[N-ii+1], N-ii+1)
        
        vals!(B[N-ii + 1], c, s)
        idx!(B[N-ii+1], N-ii+1)
    end

end
init_state{T}(state::ShiftType{T}) = QCB_factorization(state)

# If there is an issue, this function can be used to resetart the algorithm
# could be merged with init_state?
function restart{T}(state::ShiftType{T})
    # try again
    init_state(state)
    
    for i in 1:state.N
        state.REIGS[i] = state.IEIGS[i] = zero(T)
    end
    state.ctrs.zero_index = 0
    state.ctrs.start_index = 1
    state.ctrs.stop_index = state.N - 1
    state.ctrs.it_count = 0
    state.ctrs.tr = state.N - 2
end


### Related to decompostion QR into QC(B + ...)


## we need to find A[k:k+2, k:k+1] for purposes of computing eigenvalues, either
## to give the shifts or to find the roots after deflation.
##
## fill A[k:k+2, k:k+1] k in 2:N
## updates state.A
##
# We look for r_j,k. Depending on |j-k| there are different amounts of work
# we have wk = (B * D + e1 y^t) * ek = B* D * ek + e1 yk; we consider B*D* ek only B1 ... Bk D ek applies
#

## If we decompose as C'*B*D = Z, instead of D*C'*B we get
# julia> Bi*Bj*Bk* D * [0,0,1,0]
# julia> rotm(bi1, bi2, 1, 4) * rotm(bj1, bj2, 2, 4) * rotm(bk1, bk2, 3, 4) *D * [0,0,1,0]
# 4-element Array{SymPy.Sym,1}
# ⎡       ___ ___ ⎤
# ⎢bk₁⋅dk⋅bi₂⋅bj₂ ⎥
# ⎢               ⎥
# ⎢        ___ ___⎥
# ⎢-bk₁⋅dk⋅bi₁⋅bj₂⎥
# ⎢               ⎥
# ⎢         ___   ⎥
# ⎢  bk₁⋅dk⋅bj₁   ⎥
# ⎢               ⎥
# ⎣    bk₂⋅dk     ⎦
# which gives [what, wj, wk, wl]

# For rkk, we have Ck * W = [rkk, 0]
# @vars ck1 ck2 what w1
# u = rotm(ck1, ck2, 1,2) * [what, w1]
# u[1](what => solve(u[2], what)[1]) |> simplify
#     ⎛    ___      2⎞ 
# -w₁⋅⎝ck₁⋅ck₁ + ck₂ ⎠ 
# ───────────────────── this is rkk = -w1/ck_s 
#          ck₂    
# 
# For r[k-1, k] we need to do more work. We need [what_{k-1}, w_k, w_{k+1}], where w_k, w_{k+1} found from B values as above.
#
# julia> @vars ck1 ck2 cj1 cj2 what w w1
# (ck1, ck2, cj1, cj2, what, w, w1)

# julia> u = rotm(ck1, ck2, 2, 3) * rotm(cj1, cj2, 1, 3) * [what, w, w1]  # C^*_{k} * C^*{k-1} * W = [r_{k-1,k}, r_{k,k}, 0]
# 3-element Array{SymPy.Sym,1}
# ⎡        cj₁⋅ŵ - cj₂⋅w         ⎤
# ⎢                               ⎥
# ⎢                __         ⎥
# ⎢cj₂⋅ck₁⋅ŵ + ck₁⋅w⋅cj₁ - ck₂⋅w₁⎥
# ⎢                               ⎥
# ⎢               ___      ___⎥
# ⎣cj₂⋅ck₂⋅ŵ + ck₂⋅w⋅cj₁ + w₁⋅ck₁⎦



# julia> u[1](what => solve(u[3], what)[1]) |> simplify
#      2                       
#   cj₁ ⋅w   cj₁⋅ck₁⋅w₁        
# - ────── - ────────── - cj₂⋅w
#    cj₂      cj₂⋅ck₂

## or -(w + cj * ck/sk * w1) / sj

#
# For r_{k-2,k} we need to reach back one more step
# C^*_{k} * C^*{k-1} * C^*_{k-2} W = [r_{k-2,k} r_{k-1,k}, r_{k,k}, 0]
#
# julia> @vars ck1 ck2 cj1 cj2 ci1 ci2 what wm1 w w1
# julia> u = rotm(ck1, ck2, 3, 4) * rotm(cj1, cj2, 2, 4) * rotm(ci1, ci2, 1, 4) * [what, wm1, w, w1]
# julia> u[1](what => solve(u[4], what)[1]) |> simplify
#      2                                        
#   ci₁ ⋅wm₁   ci₁⋅cj₁⋅w    ci₁⋅ck₁⋅w₁          
# - ──────── - ───────── - ─────────── - ci₂⋅wm₁
#     ci₂       ci₂⋅cj₂    ci₂⋅cj₂⋅ck₂                  
#
# of -(wm1 + (ci*cj/sj)*w + (ci*ck) / (sj * sk) * w1) / si
#
# This will have problems if any of si, sj or sk are 0. This happens if the
# Ct[k] become trivial. Theorem 4.1 ensures this can't happen mathematically
# though numerically, this is a different matter. The bound involves 1/||p|| which can be smaller than machine precision for, say, Wilknson(20)
#


# D values are only for ComplexRealSingleShift
getD(state::ComplexRealSingleShift, k::Int) = state.D[k]
getD{T}(state::ShiftType{T}, k::Int) = one(T)

function diagonal_block{T}(state::ShiftType{T}, k)
    k >= 2 && k <= state.N || error("$k not in [2,n]")

    A = state.A 
    R = state.R # temp storage

    Q,Ct,B = state.Q, state.Ct, state.B
    if k == 2
        # rkk = -w_{k+1} / ck_s

        
        Bj_c, Bj_s = vals(B[k-1]); Bk_c, Bk_s = vals(B[k])
        Cj_c, Cj_s = vals(Ct[k-1]); Ck_c, Ck_s = vals(Ct[k])
        Qj_c, Qj_s = vals(Q[k-1]);  Qk_c, Qk_s = vals(Q[k])

        
        # # here we only need [r11 r12; 0 r22]
        # k=2 this is r_kk, r_k-1,k

        # for k
        wl =  Bk_s
        wk =  conj(Bj_c) * Bk_c


        R[2,2] = - wl / Ck_s
        
        # r_{k-1,k} =  -(wk + cj_c * conj(ck_c) / ck_s *wl)/cj_s
        R[1,2] = - (wk + Cj_c * conj(Ck_c) / Ck_s * wl) / Cj_s

        # for k - 1 we have (l=k)
        wl = Bj_s
        R[1,1] = - wl / Cj_s
        R[2,1] = complex(zero(T))

        alpha, beta = getD(state, k-1), getD(state, k)
        R[1,1] *= alpha; R[1,2] *= alpha
        R[2,1] *= beta; R[2,2] *= beta

# 3×2 Array{SymPy.Sym,2}
# ⎡                           ___⎤
# ⎢R₁₁⋅qjc  R₁₂⋅qjc - R₂₂⋅qkc⋅qjs⎥
# ⎢                              ⎥
# ⎢                           ___⎥
# ⎢R₁₁⋅qjs  R₁₂⋅qjs + R₂₂⋅qkc⋅qjc⎥
# ⎢                              ⎥
# ⎣   0            R₂₂⋅qks       ⎦

        A[1,1] = R[1,1] * Qj_c
        A[2,1] = R[1,1] * Qj_s
        A[1,2] = R[1,2] * Qj_c - R[2,2] * Qk_c * conj(Qj_s)
        A[2,2] = R[1,2] * Qj_s + R[2,2] * Qk_c * conj(Qj_c)

    else
        
        Bi_c, Bi_s = vals(B[k-2]);  Bj_c, Bj_s = vals(B[k-1]);  Bk_c, Bk_s = vals(B[k])
        Ci_c, Ci_s = vals(Ct[k-2]); Cj_c, Cj_s = vals(Ct[k-1]); Ck_c, Ck_s = vals(Ct[k])
        Qi_c, Qi_s = vals(Q[k-2]);  Qj_c, Qj_s = vals(Q[k-1]);  Qk_c, Qk_s = vals(Q[k])

        
        # for k
        wl =   Bk_s
        wk =  conj(Bj_c) * Bk_c
        wj = - conj(Bi_c) * conj(Bj_s) * Bk_c
        
        R[3,2] = - wl / Ck_s
        R[2,2] = - (wk + Cj_c * conj(Ck_c) / Ck_s * wl) / Cj_s
        # -(wj + ci_c * conj(cj_c) / cj_s * wk + ci_c * conj(ck_c) / (cj_s * ck_s) * wl)/ci_s
        R[1,2] = -(wj + Ci_c * conj(Cj_c) / Cj_s * wk +
                   Ci_c * conj(Ck_c) / (Cj_s * Ck_s) * wl) / Ci_s

        # downshift C indexes l->k; k->j; j->i; but keep w's (confusing)
        wl =  Bj_s
        wk =  conj(Bi_c) * Bj_c
        R[2,1] = - wl / Cj_s
        R[1,1] = - (wk + Ci_c * conj(Cj_c) / Cj_s * wl) / Ci_s
        R[3,1] = zero(T)

        alpha, beta, gamma = getD(state, k-2), getD(state, k-1), getD(state, k)
        R[1,1] *= alpha; R[1,2] *= alpha
        R[2,1] *= beta;  R[2,2] *= beta
        R[3,1] *= gamma; R[3,2] *= gamma
        
# make Qs from multiplying rotators
# make Rs = [Sym("r$i$j") for i in 1:5, j in 1:5] |> triu
# julia> (Qs * Rs)[2:4, 2:3] ## but indexing of r's is off! j-1 needed
# 3×2 Array{SymPy.Sym,2}
# ⎡                  ___                    ___               ___⎤
# ⎢q1s⋅r₁₂ + q2c⋅r₂₂⋅q1c  q1s⋅r₁₃ + q2c⋅r₂₃⋅q1c - q2s⋅q3c⋅r₃₃⋅q1c⎥
# ⎢                                                              ⎥
# ⎢                                                  ___         ⎥
# ⎢       q2s⋅r₂₂                  q2s⋅r₂₃ + q3c⋅r₃₃⋅q2c         ⎥
# ⎢                                                              ⎥
# ⎣          0                            q3s⋅r₃₃                ⎦
# 2×2 Array{SymPy.Sym,2}
# ⎡                  ___                    ___           ___ ___⎤
# ⎢R₁₂⋅qis + R₂₂⋅qjc⋅qic  R₁₃⋅qis + R₂₃⋅qjc⋅qic - R₃₃⋅qkc⋅qic⋅qjs⎥
# ⎢                                                              ⎥
# ⎢                                                  ___         ⎥
# ⎣       R₂₂⋅qjs                  R₂₃⋅qjs + R₃₃⋅qkc⋅qjc         ⎦


        A[1,1] = R[1,1] * Qi_s + R[2,1] * conj(Qi_c) * Qj_c
        A[2,1] = R[2,1] * Qj_s
        A[1,2] = R[1,2] * Qi_s + R[2,2] * conj(Qi_c) * Qj_c - R[3,2] * conj(Qi_c) * conj(Qj_s) * Qk_c
        A[2,2] = R[2,2] * Qj_s + R[3,2] * conj(Qj_c) * Qk_c
        
    end

    false 
end

##################################################

# [a11 - l a12; a21 a22] -> l^2 -2 * (tr(A)/2) l + det(A)
# so we use b = tr(A)/2 for qdrtc routing
function eigen_values{T}(state::RealDoubleShift{T})

    a11, a12 = state.A[1,1], state.A[1,2]
    a21, a22 = state.A[2,1], state.A[2,2]

    b = (a11 + a22) * (0.5)  
    c = a11 * a22 - a12 * a21
    
    state.e1[1], state.e1[2], state.e2[1], state.e2[2] = qdrtc(one(T), b, c)
    complex(state.e1[1], state.e1[2]), complex(state.e2[1], state.e2[2])            
end    

# from `modified_quadratic.f90`
function eigen_values{T}(state::ShiftType{T})

    a11, a12 = state.A[1,1], state.A[1,2]
    a21, a22 = state.A[2,1], state.A[2,2]

    tr = a11 + a22
    detm = a11 * a22 - a21 * a12
    disc = sqrt(tr * tr - 4.0 * detm)

    u = abs(tr + disc) > abs(tr - disc) ? tr + disc : tr - disc
    if iszero(u)
        state.e1[1], state.e1[2] = zero(T), zero(T)
        state.e2[1], state.e2[2] = zero(T), zero(T)
    else
        e1 = u / 2.0
        e2 = detm / e1
        state.e1[1], state.e1[2] = real(e1), imag(e1)
        state.e2[1], state.e2[2] = real(e2), imag(e2)
    end

    complex(state.e1[1], state.e1[2]), complex(state.e2[1], state.e2[2])                
end    


##################################################
## Deflation
## when a Q[k] matrix become a "D" matrix, we deflate. This is checked by the sine term being basically 0.
function check_deflation{T}(state::ShiftType{T}, tol = eps(T))
    for k in state.ctrs.stop_index:-1:state.ctrs.start_index
        if abs(vals(state.Q[k])[2]) <= tol
            deflate(state, k)
            return
        end
    end
end

# deflate a term
# turn on `show_status` to view sequence
function deflate{T}(state::ShiftType{T}, k)

    # make a D matrix
    c,s = vals(state.Q[k])
    vals!(state.Q[k], c, zero(T)) # zero out s, will renormalize c

    # shift zero counter
    state.ctrs.zero_index = k      # points to a matrix Q[k] either RealRotator(-1, 0) or RealRotator(1, 0)
    state.ctrs.start_index = k + 1

    # reset counter
    state.ctrs.it_count = 1
end


# deflate a term
# deflation for ComplexReal is different, as
# we replace Qi with I and move diagonal part into D
function deflate{T}(state::ComplexRealSingleShift{T}, k)

    # when we deflate here we want to leave Q[k] = I and
    # move Dk matrix over to merge with D
    # we do this by m
    # Qi           Qi              Dk Qi
    #   Qj     ->     Dk Qj    -->    Dk  Qj  and so on until we get to start_index
    #     Dk Ik         Dk  Ik           Dk Ik
    #
    # then the Dk's are collected into [alpa 0; I; 0 conj(alpha)] (start,k+1)

    alpha, s = vals(state.Q[k])
    vals!(state.Q[k], one(Complex{T}), zero(T)) # I

    cascade(state.Q, state.D, alpha, k, state.ctrs.stop_index) 
    
    # shift zero counter
    state.ctrs.zero_index = k      # points to a matrix Q[k] either RealRotator(-1, 0) or RealRotator(1, 0)
    state.ctrs.start_index = k + 1

    # reset counter
    state.ctrs.it_count = 1
end


##################################################



