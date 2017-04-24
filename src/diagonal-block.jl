# D values are only for ComplexRealSingleShift
#getD(state::ComplexRealSingleShift, k::Int) = state.D[k]
#getD{T}(state::ShiftType{T}, k::Int) = one(T)


### Related to decompostion QR into QC(B + ...)


## we need to find A[k:k+2, k:k+1] for purposes of computing eigenvalues, either
## to give the shifts or to find the roots after deflation.
##
## fill A[k:k+2, k:k+1] k in 2:N
## updates state.A
##
# We look for r_j,k. Depending on |j-k| there are different amounts of work
# we have wk = (B + e1 y^t) * ek = B * ek + e1 yk; we consider B* ek only B1 ... Bk D ek applies
#

# julia> Bi*Bj*Bk * [0,0,1,0]
# julia> rotm(bi1, bi2, 1, 4) * rotm(bj1, bj2, 2, 4) * rotm(bk1, bk2, 3, 4) * [0,0,1,0]
# 4-element Array{SymPy.Sym,1}
# ⎡       ___ ___ ⎤
# ⎢bk₁⋅bi₂⋅bj₂ ⎥
# ⎢               ⎥
# ⎢        ___ ___⎥
# ⎢-bk₁⋅bi₁⋅bj₂⎥
# ⎢               ⎥
# ⎢         ___   ⎥
# ⎢  bk₁⋅bj₁   ⎥
# ⎢               ⎥
# ⎣    bk₂     ⎦
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
# This will have problems if any of si, sj or sk are 0. This happens
# if the Ct[k] become trivial. Theorem 4.1 ensures this can't happen
# mathematically though numerically, this is a different matter. The
# lower bound involves 1/||p|| which can be smaller than machine
# precision for, say, Wilknson(20)
#

# R = Ct* B
# This finds R[(k-2):k,(k-1):k] from C and B,
# If D is necessary, it must be done separately
# Here k-2 <= j <= k
function Rjk(Ct, B, j, k)
    delta = k -j
    delta < 0 && return 0
    0 <= delta <= 2 || error("k-j in 0,1,2")

    if delta == 0 # (k,k) case is easiest
        # rkk = -w_{k+1} / ck_s
        wl = B[k].s
        return -wl / Ct[k].s
    elseif delta == 1
        j = k-1
        # r_{k-1,k} =  -(wk + cj_c * conj(ck_c) / ck_s *wl)/cj_s
        wl = B[k].s
        wk = conj(B[j].c) * B[k].c
        return -(wk + Ct[j].c * conj(Ct[k].c) / Ct[k].s * wl) / Ct[j].s
    else
        i,j = k-2, k-1
        wl = B[k].s
        wk = conj(B[j].c) * B[k].c
        wj = - conj(B[i].c) * conj(B[j].s) * B[k].c
        # -(wj + ci_c * conj(cj_c) / cj_s * wk + ci_c * conj(ck_c) / (cj_s * ck_s) * wl)/ci_s
        return -(wj + Ct[i].c * conj(Ct[j].c) / Ct[j].s * wk +
                 Ct[i].c * conj(Ct[k].c) / (Ct[j].s * Ct[k].s) * wl) / Ct[i].s
    end
end


# Simple type to help with getting elements of R
abstract type RType end
type RTypeDouble <: RType
    Ct
    B
end
type RTypeSingle <: RType
    Ct
    B
    D
end

Base.getindex(a::RTypeDouble, j, k) = Rjk(a.Ct, a.B, j, k)
Base.getindex(a::RTypeSingle, j, k) = a.D[j] * Rjk(a.Ct, a.B, j, k)

function diagonal_block{T, St, Tw}(state::FactorizationType{T, St, Val{:NoPencil}, Tw}, k)
    k >= 2 && k <= state.N || error("$k not in [2,n]")

    R = state.R # temp storage

    Q,Ct,B = state.Q, state.Ct, state.B
    
    V = St == Val{:SingleShift} ? RTypeSingle(state.Ct, state.B, state.D) : RTypeDouble(state.Ct, state.B)
    
    if k == 2
        for i in 1:2, j in 1:2
            R[i,j] = V[i,j]
        end
    else
        for i in 1:3, j in 1:2
            R[i,j] = V[k - 3 + i, k - 2 + j]
        end
    end

    compute_QR(state, R, k)

    false 
end
## Has pencil!!

# julia> s = VWi[2:4, 3:4]
# 3×2 Array{SymPy.Sym,2}:
#  v23/W33 - W23*v22/(W22*W33)  v24/W44 - W34*v23/(W33*W44) + v22*(W23*W34/(W33*W44) - W24/W44)/W22
#                      v33/W33                                          v34/W44 - W34*v33/(W33*W44)
#                            0                                                              v44/W44


function diagonal_block{T, St, Tw}(state::FactorizationType{T, St, Val{:HasPencil}, Tw}, k)

    k >= 2 && k <= state.N || error("$k not in [2,n]")

    R = state.R # temp storage

    W = St == Val{:SingleShift} ? RTypeSingle(state.Ct1, state.B1, state.D1) : RTypeDouble(state.Ct1, state.B1)    
    V = St == Val{:SingleShift} ? RTypeSingle(state.Ct, state.B, state.D) : RTypeDouble(state.Ct, state.B)



    if k == 2
        j = k - 1
        R[1,1] = V[j,j] / W[j,j]
        R[1,2] = V[j,k]/W[k,k] - W[j,k]*V[j,j] / (W[j,j]*W[k,k])
        R[2,1] = zero(T)
        R[2,2] = V[k,k] / W[k,k]

    else
        
        i,j = k-2, k-1
        R[1,1] = V[i,j]/W[j,j] - W[i,j] * V[i,i] / (W[i,i]*W[j,j])
        R[1,2] = V[i,k]/W[k,k] - W[j,k] * V[i,j] / (W[j,j]*W[k,k]) + V[i,i] *(W[i,j]*W[j,k]/(W[i,i]*W[j,j]) - W[i,k]/W[k,k])/W[i,i]
        R[2,1] = V[j,j] / W[j,j]
        R[2,2] = V[j,k]/W[k,k] - W[j,k]*V[j,j] - W[j,k]*V[j,j] / (W[i,i]*W[j,j])
        R[3,1] = zero(T)
        R[3,2] = V[k,k] / W[k,k]
    end        

    # rotate by Q
    compute_QR(state, R, k)
        

    false 
end



# allows treatment of complex and real case at same time
getD{T, P, Tw}(state::FactorizationType{T, Val{:SingleShift}, P, Tw}, k) = state.D[k]
getD{T, P, Tw}(state::FactorizationType{T, Val{:DoubleShift}, P, Tw}, k) = one(T)

## we compute R from Ct, B
## we compute QR from R and Q [D],
## This allows us to handle twisting separately
## Pencil doesn't effect this code, that is the R part
function compute_QR{T,St, P}(state::FactorizationType{T, St, P, Val{:NotTwisted}}, R, k)
    A = state.A
    Q = state.Q

    if k == 2
       
# 3×2 Array{SymPy.Sym,2}
# ⎡                           ___⎤
# ⎢R₁₁⋅qjc  R₁₂⋅qjc - R₂₂⋅qkc⋅qjs⎥
# ⎢                              ⎥
# ⎢                           ___⎥
# ⎢R₁₁⋅qjs  R₁₂⋅qjs + R₂₂⋅qkc⋅qjc⎥
# ⎢                              ⎥
# ⎣   0            R₂₂⋅qks       ⎦

        ## combine with Qj_c
        Qj_c, Qj_s = vals(Q[k-1]);  Qk_c, Qk_s = vals(Q[k])

        A[1,1] = R[1,1] * Qj_c
        A[2,1] = R[1,1] * Qj_s
        A[1,2] = R[1,2] * Qj_c - R[2,2] * Qk_c * conj(Qj_s)
        A[2,2] = R[1,2] * Qj_s + R[2,2] * Qk_c * conj(Qj_c)
    else
     


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
        Qi_c, Qi_s = vals(Q[k-2]);  Qj_c, Qj_s = vals(Q[k-1]);  Qk_c, Qk_s = vals(Q[k])        
        
        A[1,1] = R[1,1] * Qi_s + R[2,1] * conj(Qi_c) * Qj_c
        A[2,1] = R[2,1] * Qj_s
        A[1,2] = R[1,2] * Qi_s + R[2,2] * conj(Qi_c) * Qj_c - R[3,2] * conj(Qi_c) * conj(Qj_s) * Qk_c
        A[2,2] = R[2,2] * Qj_s + R[3,2] * conj(Qj_c) * Qk_c

    end
end



##################################################

# [a11 - l a12; a21 a22] -> l^2 -2 * (tr(A)/2) l + det(A)
# so we use b = tr(A)/2 for qdrtc routing
function eigen_values{T,P, Tw}(state::FactorizationType{T,Val{:DoubleShift}, P, Tw})

    a11, a12 = state.A[1,1], state.A[1,2]
    a21, a22 = state.A[2,1], state.A[2,2]

    b = (a11 + a22) * (0.5)  
    c = a11 * a22 - a12 * a21
    
    state.e1[1], state.e1[2], state.e2[1], state.e2[2] = qdrtc(one(T), b, c)
    complex(state.e1[1], state.e1[2]), complex(state.e2[1], state.e2[2])            
end    

# from `modified_quadratic.f90`
function eigen_values{T,P, Tw}(state::FactorizationType{T,Val{:SingleShift}, P, Tw})

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

