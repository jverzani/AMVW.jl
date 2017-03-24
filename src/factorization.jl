## 
## initial factorization

# ## This is working for RDS
# function init_state_XXX{T}(state::ShiftType{T})

#     N, ps= state.N, state.POLY
#     par = iseven(N) ? one(T) : -one(T)
    
#     Q, Ct, B = state.Q, state.Ct, state.B

     
#     for ii = 1:(N-1)
#         vals!(Q[ii], zero(T), one(T))
#         idx!(Q[ii], ii)
#     end
#     vals!(Q[N], one(T), zero(T))
#     idx!(Q[N], N)



#     # ## Working, but not quite what is in DFCC code (there -par*ps[N], par*one(T); C is -conj(c), -s
#     # c, s, temp = givensrot(par * ps[N], -par * one(T))
#     # vals!(Ct[N], conj(c), -s); idx!(Ct[N], N)


#     c, s, temp = givensrot(par * ps[N], -par * one(T))
#     vals!(Ct[N], conj(c), -s); idx!(Ct[N], N)

#     vals!(B[N], -par * conj(s), par * c)
#     idx!(B[N], N)
    
  
        
    
#     for ii in 2:N
#         c, s, temp = givensrot(-ps[ii-1], temp)
#         vals!(Ct[N-ii + 1], conj(c), -s)
#         idx!(Ct[N-ii+1], N-ii+1)
        
#         vals!(B[N-ii + 1], c, s)
#         idx!(B[N-ii+1], N-ii+1)
#     end

# end


## This is working for RDS
function init_state{T}(state::ShiftType{T})

    N, ps= state.N, state.POLY
    par = iseven(N) ? one(T) : -one(T)
    
    Q, Ct, B = state.Q, state.Ct, state.B

    for ii = 1:(N-1)
        if isa(state, RealDoubleShift)
            vals!(Q[ii], zero(T), one(T))
        else
            vals!(Q[ii], zero(T), zero(T), one(T))
        end
        idx!(Q[ii], ii)
    end
    if isa(state, RealDoubleShift)
        vals!(Q[N], one(T), zero(T))
    else
        vals!(Q[N], one(T), zero(T), zero(T))
    end
    idx!(Q[N], N)


    ## Working, but not quite what is in DFCC code
    ## there -par*ps[N], par*one(T); C is -conj(c), -s
    ## B[N] = par, -par...
    ## Here we use: C' * B * D = Z (and not D C' B = Z)
    c, s, temp = givensrot(par * ps[N], -par * one(T))
    vals!(Ct[N], conj(c), -s); idx!(Ct[N], N)

    if isa(c, Real)
        vals!(B[N], -par * conj(s), par * c) 
        idx!(B[N], N)
    else
        if isreal(c)
            r = real(c)
            phi = complex(one(T), zero(T))
        else
            r = norm(c)
            phi = c / r
        end
        
        vals!(B[N], -par * s * phi, par * r)
        idx!(B[N], N)
        state.D[N] = conj(phi)
        state.D[N+1] = phi
    end

    
    
    for ii in 2:N
        c, s, temp = givensrot(-ps[ii-1], temp)
        vals!(Ct[N-ii + 1], conj(c), -s)
        idx!(Ct[N-ii+1], N-ii+1)
        
        vals!(B[N-ii + 1], c, s)
        idx!(B[N-ii+1], N-ii+1)
    end

end


# If there is an issue, this function can be used to resetart the algorithm
# could be merged with init_state?
function restart{T}(state::RealDoubleShift{T})
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

# ## This is working for RDS
# function init_state_WORKING{T}(state::ShiftType{T})

#     N, ps= state.N, state.POLY
#     par = iseven(N) ? one(T) : -one(T)
    
#     Q, Ct, B = state.Q, state.Ct, state.B
    
#     for ii = 1:(N-1)
#         vals!(Q[ii], zero(T), one(T))
#         idx!(Q[ii], ii)
#     end
#     vals!(Q[N], one(T), zero(T))
#     idx!(Q[N], N)


#     ## Working, but not quite what is in DFCC code (there -par*ps[N], par*one(T); C is -conj(c), -s
#     c, s, temp = givensrot(-ps[N], -one(T))
#     vals!(Ct[N], -par*conj(c), -par*s); idx!(Ct[N], N)
    
#     if isa(c, Real)
#         vals!(B[N], -conj(s), -c) 
#         idx!(B[N], N)
#     else
#         if isreal(c)
#             r = real(c)
#             phi = complex(one(T), zero(T))
#         else
#             r = norm(c)
#             phi = conj(c) / r
#         end

#         vals!(B[N], -s * conj(phi), -r)
#         idx!(B[N], N)
#         state.D[N] = conj(phi)
#         state.D[N+1] = phi
#     end

    
    
#     for ii in 2:N
#         c, s, temp = givensrot(-ps[ii-1], temp)
#         vals!(Ct[N-ii + 1], conj(c), -s)
#         idx!(Ct[N-ii+1], N-ii+1)
        
#         vals!(B[N-ii + 1], c, s)
#         idx!(B[N-ii+1], N-ii+1)
#     end

# end


# ## Do I need to changes for CSS????
# function init_state{T}(state::ShiftType{T})

#     N, ps= state.N, state.POLY
#     par = iseven(N) ? one(T) : -one(T)
    
#     Q, Ct, B = state.Q, state.Ct, state.B
    
#     for ii = 1:(N-1)
#         vals!(Q[ii], zero(T), one(T))
#         idx!(Q[ii], ii)
#     end
#     vals!(Q[N], one(T), zero(T))
#     idx!(Q[N], N)


#     ## Follow factor.f90
#     c, s, temp = givensrot(par * ps[N], -par * one(T))
#     vals!(Ct[N], -conj(c), -s); idx!(Ct[N], N)
    
#     if isa(c, Real)
#         vals!(B[N], -conj(s), -c) 
#         idx!(B[N], N)
#     else
#         if isreal(c)
#             r = real(c)
#             phi = complex(one(T), zero(T))
#         else
#             r = norm(c)
#             phi = conj(c) / r
#         end

#         vals!(B[N], -s * conj(phi), -r)
#         idx!(B[N], N)
#         state.D[N] = conj(phi)
#         state.D[N+1] = phi
#     end

    
    
#     for ii in 2:N
#         c, s, temp = givensrot(-ps[ii-1], temp)
#         vals!(Ct[N-ii + 1], conj(c), -s)
#         idx!(Ct[N-ii+1], N-ii+1)
        
#         vals!(B[N-ii + 1], c, s)
#         idx!(B[N-ii+1], N-ii+1)
#     end

# end



# function init_stateA{T}(state::ComplexSingleShift{T})
#     N, ps= state.N, state.POLY
#     p = iseven(N) ? one(T) : -one(T)

#     Q, Ct, B = state.Q, state.Ct, state.B
    
#     for ii = 1:(N-1)
#         vals!(Q[ii], zero(T), zero(T), one(T)); idx!(Q[ii], ii)
#     end

#     vals!(Q[N], one(T), zero(T), zero(T)); idx!(Q[N], N)

#     # # play with signs here.
#     # s = iseven(N) ? one(T) : -one(T)

#     # a, b, temp = givensrot(-ps[N], -one(T))
#     # println("givens rot $(-ps[N]), -1 -> $a, $b, $temp")


#     # vals!(Ct[N], -a, -b); idx!(Ct[N], N)

#     # abarhat = conj(a)/norm(a)  ## udpate D!
#     # state.D[N] = abarhat
    
#     # vals!(B[N], -s * b * abarhat, -s * norm(a)) 
#     # idx!(B[N], N)


    
#         c, s, temp = givensrot(-1 * ps[N], 1 *one(T))
#     #    c, s, temp = givensrot( p * ps[N], -p  * one(T))
# #    c, s, temp = givensrot( p * ps[N], - p * one(T))        
#     # println("givens rot $(p*ps[N]), $(-p) -> $c, $s, $temp")


# #    vals!(Ct[N], conj(c), -s); idx!(Ct[N], N) #  C is (c, s)
#     vals!(Ct[N], -conj(c), -s); idx!(Ct[N], N) #  C is (c, s)    

#     # B would be C(s, c) but, need to rotate so that complex value is the cosine
#     if iszero(imag(c))
#         r = real(c)
#         phi = complex(one(T), zero(T))
#     else
#         r = norm(c)
#         phi = conj(c) / r  ##
#     end
    
# #    vals!(B[N],   -s *  conj(phi),  r)
#     vals!(B[N],   -s *  conj(phi),  -r)         
#     idx!(B[N], N)

#     # D is inverse of [phi 0; 0 conj(phi)] -> [conj(phi) 0; 0 phi]
#     state.D[N] = conj(phi)
#     state.D[N+1] =  phi
# #    state.D[N] =  -p*conj(phi)
# #    state.D[N+1] =  -p*phi


#     # check B = C D' Z (We store D'
#     # Z = [0 -1; 1 0]; D = [state.D[N] 0; 0 state.D[N+1]]'; C = [c -s; s c]
#     # Bc = [-s*conj(phi) -r; r -s*phi]
#     # println("B,CDZ")
#     # println(Bc)
#     # println(C*D*Z)

    
#     for ii in 2:N
#         temp1 = temp
#         c, s, temp = givensrot(-ps[ii-1], temp)
# #        println("XXX givens rot $(-ps[ii-1]), $temp1 -> $c, $s, $temp")        
#         vals!(Ct[N-ii + 1], conj(c), -s); idx!(Ct[N-ii+1], N-ii+1)
#         vals!(B[N-ii + 1], c, s); idx!(B[N-ii+1], N-ii+1)
#     end

# end

##
### Related to decompostion QR into QC(B + ...)


## we need to find A[k:k+2, k:k+1] for purposes of computing eigenvalues, either
## to give the shifts or to find the roots after deflation.
##
## fill A[k:k+2, k:k+1] k in 2:N
## updates state.A
##
# We look for r_j,k. Depending on |j-k| there are different amounts of work
# we have wk = (B + e1 y^t) * ek = B*ek + e1 yk; we consider B * ek only B1 ... Bk ek applies
#
# julia> @vars bk1 bk2 bj1 bj2 bi1 bi2
# julia> rotm(bi1, bi2, 1, 4) * rotm(bj1, bj2, 2, 4) * rotm(bk1, bk2, 3, 4) * [0, 0, 1, 0]  # B_{k-2} * B_{k-1} * B_k * ek = W
# 4-element Array{SymPy.Sym,1}
# ⎡bi₂⋅bj₂⋅bk₁ ⎤
# ⎢            ⎥
# ⎢         ___⎥
# ⎢-bj₂⋅bk₁⋅bi₁⎥
# ⎢            ⎥
# ⎢      ___   ⎥
# ⎢  bk₁⋅bj₁   ⎥
# ⎢            ⎥
# ⎣    bk₂     ⎦
# which gives W = [what_{k-2} w_{k-1} w_k w_{k+1}]

## If we decompose as C'*D * B = Z, instead of D*C'*B we get
# julia> D*Bi*Bj*Bk*[0,0,1,0]
# 4-element Array{SymPy.Sym,1}
# ⎡       ___ ___ ⎤
# ⎢bk₁⋅di⋅bi₂⋅bj₂ ⎥
# ⎢               ⎥
# ⎢        ___ ___⎥
# ⎢-bk₁⋅dj⋅bi₁⋅bj₂⎥
# ⎢               ⎥
# ⎢         ___   ⎥
# ⎢  bk₁⋅dk⋅bj₁   ⎥
# ⎢               ⎥
# ⎣    bk₂⋅dl     ⎦

# For rkk, we have Ck * W = [rkk, 0]
# @vars ck1 ck2 what w1
# u = rotm(ck1, ck2, 1,2) * [what, w1]
# u[1](what => solve(u[2], what)[1]) |> simplify
#     ⎛    ___      2⎞ 
# -w₁⋅⎝ck₁⋅ck₁ + ck₂ ⎠ 
# ─────────────────────
#          ck₂    

#    ⎛  2     2⎞ 
# -w₁⋅⎝c₁  + c₂ ⎠ 
# ────────────────  # this is rkk = -w1/c2 = -bk2/ck2
#        c₂ 

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
# function diagonal_blockA{T}(state::RealDoubleShift{T}, k)
#     k >= 2 && k <= state.N || error("$k not in [2,n]")

#     A = state.A
#     R = state.R

#     if k == 2
#         # # here we only need [r11 r12; 0 r22], so only use top part of R
#         for j in 1:2
#             ck1, ck2 = vals(state.Ct[k - (2-j)])
#             w1 = vals(state.B[k - (2-j)])[2]
#             R[j,j] = - w1 / ck2
#         end
        
#         cj1, cj2 = vals(state.Ct[k-1])
#         ck1, ck2 = vals(state.Ct[k])

#         ## should always be non zero, but...
#         (iszero(cj2)  || iszero(ck2)) && return false
        
#         w = vals(state.B[k-1])[1] * vals(state.B[k])[1]
#         w1 = vals(state.B[k])[2]
#         val = -(1/cj2) * (w + (cj1 * ck1 / ck2) * w1)
#         R[1,2] = val

#         q11, q12 = vals(state.Q[k-1]); q21, q22 = vals(state.Q[k])
#         A[1,1] = q11 * R[1,1]
#         A[1,2] = q11 * R[1,2] - q12 * q21 * R[2,2]
#         A[2,1] = q12 * R[1,1]
#         A[2,2] = q11 * q21 * R[2,2] + q12 * R[1,2]

        
                   
#     else  

#         ci, si = vals(state.Ct[k-2])

#         ## shouldn't happen, but can
#         iszero(si) && return false

        
#         Qk_2, Qk_1, Qk = vals(state.Q[k-2]), vals(state.Q[k-1]), vals(state.Q[k])
#         wk_1, wk, wk1 = vals(state.B[k-2])[2], vals(state.B[k-1])[2], vals(state.B[k])[2]

#         # r_kk
#         for j in 1:2
#             K = k - (2-j)
#             ck, sk = vals(state.Ct[K])
#             w1 = vals(state.B[K])[2]
#             R[j+1,j] = - w1 / sk
#         end

#         #R_{k,k-1}
#         for j in 1:2
#             K = k - (2-j)
#             cj, sj = vals(state.Ct[K-1])
#             ck, sk = vals(state.Ct[K])
#             w = vals(state.B[K-1])[1] * vals(state.B[K])[1]
#             w1 = vals(state.B[K])[2]
#             # need ck2, cj2 non zero
#             val = -(w  +    (cj * ck / sk) * w1) / sj
#             R[j,j] = val
#         end

#         # R_{k, k-2}  is R[1,2]
#         wm1 = -vals(state.B[k-2])[1] * vals(state.B[k-1])[2] * vals(state.B[k])[1]
#         w = vals(state.B[k-1])[1] * vals(state.B[k])[1]
#         w1 = vals(state.B[k])[2]

#         ci, si = vals(state.Ct[k-2])
#         cj, sj = vals(state.Ct[k-1])
#         ck, sk = vals(state.Ct[k])

#         val = -(wm1 + (ci * cj /sj) * w + (ci * ck / sj / sk) * w1) / si
        
        
#         R[1,2] = val

#         #-(ci1^2 * wm1/ci2) - (ci1 * cj1 * w) / (ci2 * cj2) - (ci1 * ck1 * w1) / (ci2 * cj2 * ck2) - ci2*wm1

#         # A is Q*R, but not all Q contribute
#         # This is for k = 5
#         # julia> A[4,4]
#         # q₃ ₁⋅q₄ ₁⋅r₄₄ + q₃ ₂⋅r₃₄
        
#         # julia> A[4,5]
#         # q₃ ₁⋅q₄ ₁⋅r₄₅ - q₃ ₁⋅q₄ ₂⋅q₅ ₁⋅r₅₅ + q₃ ₂⋅r₃₅
        
#         # julia> A[5,4]
#         # q₄ ₂⋅r₄₄
        
#         # julia> A[5,5]
#         # q₄ ₁⋅q₅ ₁⋅r₅₅ + q₄ ₂⋅r₄₅

#         println(R)
        
#         A[1,1] = Qk_2[1] * Qk_1[1] * R[2,1] + Qk_2[2] * R[1,1]
#         A[1,2] = Qk_2[1] * Qk_1[1] * R[2,2] - Qk_2[1] * Qk_1[2] * Qk[1] * R[3,2] + Qk_2[2] * R[1,2]
#         A[2,1] = Qk_1[2] * R[2,1]
#         A[2,2] = Qk_1[1] * Qk[1] * R[3,2] + Qk_1[2] * R[2,2]

        
#     end
#     return true
# end

##
## rkk = - bk_s/ck_s
## rj,k =
#             ____ ____                              ____
#   bk_c⋅cj_c⋅bj_c⋅cj_c             ____   bk_s⋅cj_c⋅ck_c
# - ─────────────────── - bk_c⋅cj_s⋅bj_c - ──────────────
#           cj_s                             cj_s⋅ck_s   
#
#
# ri,k =

# using SymPy
# function rotm(c,s, i, N)
#     r = eye(Sym, N)
#     r[i:i+1, i:i+1] = [c -s; s conj(c)]
#     r
# end

# ## We solve wk = (B + e1 y^t) * ek = B*ek + e1 yk. Consdier B*ek

# @vars bi_c bi_s bj_c bj_s bk_c bk_s
# @vars what wj wk wl  # wl is w[k+1]
# @vars ci_c ci_s cj_c cj_s ck_c ck_s


# W = rotm(bi_c, bi_s, 1, 4) * rotm(bj_c, bj_s, 2, 4) * rotm(bk_c, bk_s, 3, 4) * [0, 0, 1, 0] # [what, wj, wk, wl]
# what, wj, wk, wl = W

# ## We have Ck * W = [rkk, 0] so
# u = rotm(ck_c, ck_s, 1, 2) * [what, wl]
# rkk = u[1](what => solve(u[end], what)[1]) |> simplify

# #      ⎛     ____       2⎞ 
# # -bk_s⋅⎝ck_c⋅ck_c + ck_s ⎠ 
# # ────────────────────────── = - bk_s/ck_s  (cos^2 + sin^2 = 1)
# #            ck_s

# W = [what, wk, wl]
# u = rotm(ck_c, ck_s, 2, 3) * rotm(cj_c, cj_s, 1, 3) * W  # C^*_{k} * C^*{k-1} * W  = [r_k-1,k, r_kk, 0
# rjk = u[1](what => solve(u[end], what)[1]) |> simplify

# # rjk = 
# #             ____ ____                            ____
# #   bk_c⋅cj_c⋅bj_c⋅cj_c             ____   bk_s⋅cj_c⋅ck_c
# # - ─────────────────── - bk_c⋅cj_s⋅bj_c - ──────────────
# #           cj_s                             cj_s⋅ck_s   


# W = [what, wj, wk, wl]
# u = rotm(ck_c, ck_s, 3, 4) * rotm(cj_c, cj_s, 2, 4) * rotm(ci_c, ci_s, 1, 4) * W  # C^*_{k} * C^*{k-1} * W
# rjk = u[1](what => solve(u[end], what)[1]) |> simplify

# # rik = 
# #                ____ ____                              ____ ____                ____
# # bj_s⋅bk_c⋅ci_c⋅bi_c⋅ci_c                  ____   bk_c⋅ci_c⋅bj_c⋅cj_c      bk_s⋅ci_c⋅ck_c
# # ──────────────────────── + bj_s⋅bk_c⋅ci_s⋅bi_c - ─────────────────── - ──────────────
# #           ci_s                                        ci_s⋅cj_s        ci_s⋅cj_s⋅ck_s

#function diagonal_block{T}(state::ComplexSingleShift{T}, k)
D_adjust(state::RealDoubleShift, R, k) = nothing
function D_adjust(state::ComplexSingleShift, R, k)
    println("D adjust $k")
    D = state.D
    if k == 2
        R[1,:] = state.D[k-1] * R[1,:]
        R[2,:] = state.D[k-0] * R[2,:]
    else
        println("testing 123")
        d0, d1, d2, d3 =  D[k-2], D[k-1], D[k-0], D[k+1]
        println((d0, d1, d2, d3))
        println(R)
        println("---")

        
        R[1,:] = state.D[k-2] * R[1,:]
        R[2,:] = state.D[k-1] * R[2,:]
        R[3,:] = state.D[k-0] * R[3,:]
    end
end

get_d(state::ComplexSingleShift,k) = state.D[k]
get_d{T}(state::RealDoubleShift{T},k) = one(T)

function diagonal_block{T}(state::ShiftType{T}, k)
    k >= 2 && k <= state.N || error("$k not in [2,n]")

    A = state.A # holds matrix in 1:2, 1:2
    R = state.R # temp storate

    Q,Ct,B = state.Q, state.Ct, state.B
    if k == 2
        # rkk = -w_{k+1} / ck_s


        
        Bj_c, Bj_s = vals(B[k-1]); Bk_c, Bk_s = vals(B[k])
        Cj_c, Cj_s = vals(Ct[k-1]); Ck_c, Ck_s = vals(Ct[k])
        Qj_c, Qj_s = vals(Q[k-1]);  Qk_c, Qk_s = vals(Q[k])

        
        # # here we only need [r11 r12; 0 r22], so only use top part of R
        # follow fortran code
        # k=2 this is r_kk, r_k-1,k

        # for k
        wl = get_d(state,k) * Bk_s
        wk = get_d(state, k) * conj(Bj_c) * Bk_c
        R[2,2] = - wl / Ck_s

##                             ___            
##         ___   cj₁⋅w₁⋅ck₁         ___
## - cj₁⋅w⋅cj₁ - ────────── - cj₂⋅w⋅cj₂   = 
##                  ck₂                
        
        # r_{k-1,k} =  -(wk + cj_c * conj(ck_c) / ck_s *wl)/cj_s
        R[1,2] = - (wk + Cj_c * conj(Ck_c) / Ck_s * wl) / Cj_s

        # for k - 1 we have (l=k)
        wl = get_d(state, k-1) * Bj_s
        R[1,1] = - wl / Cj_s
        R[2,1] = complex(zero(T))


# julia> rotm(qjc, qjs, 1, 3) * rotm(qkc, qks, 2, 3) * [r11 r12; 0 r22; 0 0]
# 3×2 Array{SymPy.Sym,2}
# ⎡qjc⋅r₁₁  qjc⋅r₁₂ - qjs⋅qkc⋅r₂₂⎤
# ⎢                              ⎥
# ⎢                           ___⎥
# ⎢qjs⋅r₁₁  qjs⋅r₁₂ + qkc⋅r₂₂⋅qjc⎥
# ⎢                              ⎥
# ⎣   0            qks⋅r₂₂       ⎦        

        A[1,1] = Qj_c * R[1,1]
        A[1,2] = Qj_c * R[1,2] - Qj_s * Qk_c * R[2,2]
        A[2,1] = Qj_s * R[1,1]
        A[2,2] = Qj_s * R[1,2] + conj(Qj_c) * Qk_c * R[2,2] 


    else
        
        Bi_c, Bi_s = vals(B[k-2]); Bj_c, Bj_s = vals(B[k-1]); Bk_c, Bk_s = vals(B[k])
        Ci_c, Ci_s = vals(Ct[k-2]); Cj_c, Cj_s = vals(Ct[k-1]); Ck_c, Ck_s = vals(Ct[k])
        Qi_c, Qi_s = vals(Q[k-2]); Qj_c, Qj_s = vals(Q[k-1]); Qk_c, Qk_s = vals(Q[k])

        
        # for k
        wl = get_d(state, k) * Bk_s
        wk = get_d(state, k) * conj(Bj_c) * Bk_c
        wj = - get_d(state, k) * conj(Bi_c) * conj(Bj_s) * Bk_c
        
        R[3,2] = - wl / Ck_s
        R[2,2] = - (wk + Cj_c * conj(Ck_c) / Ck_s * wl) / Cj_s
        # -(wj + ci_c * conj(cj_c) / cj_s * wk + ci_c * conj(ck_c) / (cj_s * ck_s) * wl)/ci_s
        R[1,2] = -(wj + Ci_c * conj(Cj_c) / Cj_s * wk +
                   Ci_c * conj(Ck_c) / (Cj_s * Ck_s) * wl) / Ci_s

        # downshift C indexes l->k; k->j; j->i; but keep w's (confusing)
        wl = get_d(state, k-1) * Bj_s
        wk = get_d(state, k-1) * conj(Bi_c) * Bj_c
        R[2,1] = - wl / Cj_s
        R[1,1] = - (wk + Ci_c * conj(Cj_c) / Cj_s * wl) / Ci_s
        R[3,1] = zero(T)

        ## XXX These 4 work for real, but not XXX
#        A[1,1] = Qi_c * R[1,1] - Qi_s * Qj_c * R[2,1]
#        A[1,2] = Qi_s * R[1,2] + conj(Qi_c) * Qj_c * R[2,2] - Qi_c * Qj_s * Qk_c * R[3,2]  
#        A[2,1] = Qj_s * R[2,1]
#        A[2,2] = Qj_c * Qk_c * R[3,2] + Qj_s * R[2,2]

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
        
        A[1,2] = Qi_s * R[1,1] + conj(Qi_c) * Qj_c * R[2,1]
        A[1,2] = Qi_s * R[1,2] + conj(Qi_c) * Qj_c * R[2,2] - conj(Qi_c) * Qj_s * Qk_c * R[3,2]
        A[2,1] = Qj_s * R[2,1]
        A[2,2] = Qj_s * R[2,2] + Qk_c * conj(Qj_c) * R[3,2]

    end

    false ## XXX work on flag
end

    
# function diagonal_block_WORKING{T}(state::ShiftType{T}, k)
#     k >= 2 && k <= state.N || error("$k not in [2,n]")

#     A = state.A # holds matrix in 1:2, 1:2
#     R = state.R # temp storate

#     Q,Ct,B = state.Q, state.Ct, state.B
#     if k == 2
#         Bj_c, Bj_s = vals(B[k-1]); Bk_c, Bk_s = vals(B[k])
#         Cj_c, Cj_s = vals(Ct[k-1]); Ck_c, Ck_s = vals(Ct[k])
#         Qj_c, Qj_s = vals(Q[k-1]);  Qk_c, Qk_s = vals(Q[k])

        
#         # # here we only need [r11 r12; 0 r22], so only use top part of R
#         # follow fortran code
#         R[1,1] = - Bj_s / Cj_s
#         R[2,2] = - Bk_s / Ck_s
#         R[1,2] = - (Bk_c * conj(Bj_c) + Bk_s * Cj_c * conj(Ck_c)/Ck_s) / Cj_s
#         R[2,1] = complex(zero(T))

# ## XX        D_adjust(state, R, k)

           

# # julia> rotm(qjc, qjs, 1, 3) * rotm(qkc, qks, 2, 3) * [r11 r12; 0 r22; 0 0]
# # 3×2 Array{SymPy.Sym,2}
# # ⎡qjc⋅r₁₁  qjc⋅r₁₂ - qjs⋅qkc⋅r₂₂⎤
# # ⎢                              ⎥
# # ⎢                           ___⎥
# # ⎢qjs⋅r₁₁  qjs⋅r₁₂ + qkc⋅r₂₂⋅qjc⎥
# # ⎢                              ⎥
# # ⎣   0            qks⋅r₂₂       ⎦        

#         A[1,1] = Qj_c * R[1,1]
#         A[1,2] = Qj_c * R[1,2] - Qj_s * Qk_c * R[2,2]
#         A[2,1] = Qj_s * R[1,1]
#         A[2,2] = Qj_s * R[1,2] + conj(Qj_c) * Qk_c * R[2,2] 


#     else
        
#         Bi_c, Bi_s = vals(B[k-2]); Bj_c, Bj_s = vals(B[k-1]); Bk_c, Bk_s = vals(B[k])
#         Ci_c, Ci_s = vals(Ct[k-2]); Cj_c, Cj_s = vals(Ct[k-1]); Ck_c, Ck_s = vals(Ct[k])
#         Qi_c, Qi_s = vals(Q[k-2]); Qj_c, Qj_s = vals(Q[k-1]); Qk_c, Qk_s = vals(Q[k])
        
        

        
#         R[2,1] = - Bj_s / Cj_s
#         R[3,2] = - Bk_s / Ck_s

#         R[1,1] = (-conj(Bi_c) * Bj_c - Bj_s * Ci_c * conj(Cj_c)/Cj_s) / Ci_s
#         R[2,2] = (-conj(Bj_c) * Bk_c - Bk_s * Cj_c * conj(Ck_c)/Ck_s) / Cj_s
        
        
#         R[1,2] = (conj(Bi_c)  * Bj_s * Bk_c  -
#                   conj(Bj_c) * Bk_c * Ci_c * conj(Cj_c) / Cj_s -
#                   Bk_s * Ci_c * conj(Ck_c) / (Cj_s * Ck_s)) / Ci_s
#         R[3,1] = zero(T)

#         D_adjust(state, R, k)

#         ## XXX These 4 work for real, but not XXX
# #        A[1,1] = Qi_c * R[1,1] - Qi_s * Qj_c * R[2,1]
# #        A[1,2] = Qi_s * R[1,2] + conj(Qi_c) * Qj_c * R[2,2] - Qi_c * Qj_s * Qk_c * R[3,2]  
# #        A[2,1] = Qj_s * R[2,1]
# #        A[2,2] = Qj_c * Qk_c * R[3,2] + Qj_s * R[2,2]

# # make Qs from multiplying rotators
# # make Rs = [Sym("r$i$j") for i in 1:5, j in 1:5] |> triu
# # julia> (Qs * Rs)[2:4, 2:3] ## but indexing of r's is off! j-1 needed
# # 3×2 Array{SymPy.Sym,2}
# # ⎡                  ___                    ___               ___⎤
# # ⎢q1s⋅r₁₂ + q2c⋅r₂₂⋅q1c  q1s⋅r₁₃ + q2c⋅r₂₃⋅q1c - q2s⋅q3c⋅r₃₃⋅q1c⎥
# # ⎢                                                              ⎥
# # ⎢                                                  ___         ⎥
# # ⎢       q2s⋅r₂₂                  q2s⋅r₂₃ + q3c⋅r₃₃⋅q2c         ⎥
# # ⎢                                                              ⎥
# # ⎣          0                            q3s⋅r₃₃                ⎦
        
#         A[1,2] = Qi_s * R[1,1] + conj(Qi_c) * Qj_c * R[2,1]
#         A[1,2] = Qi_s * R[1,2] + conj(Qi_c) * Qj_c * R[2,2] - conj(Qi_c) * Qj_s * Qk_c * R[3,2]
#         A[2,1] = Qj_s * R[2,1]
#         A[2,2] = Qj_s * R[2,2] + Qk_c * conj(Qj_c) * R[3,2]

#     end

#     false ## XXX work on flag
# end


# [a11 - l a12; a21 a22] -> l^2 -2 * (tr(A)/2) l + det(A)
# so we use b = tr(A)/2 for qdrtc routing
function eigen_values{T}(state::RealDoubleShift{T})

    a11, a12 = state.A[1,1], state.A[1,2]
    a21, a22 = state.A[2,1], state.A[2,2]

    b = (a11 + a22) * (0.5)  
    c = a11 * a22 - a12 * a21
    
    state.e1[1], state.e1[2], state.e2[1], state.e2[2] = qdrtc(one(T), b, c)
end    

function eigen_values{T}(state::ComplexSingleShift{T})

    a11, a12 = state.A[1,1], state.A[1,2]
    a21, a22 = state.A[2,1], state.A[2,2]
    ## XXX work this out...XXX
    e1, e2 = eigvals(state.A)
    state.e1[1], state.e1[2] = real(e1), imag(e1)
    state.e2[1], state.e2[2] = real(e2), imag(e2)    
end    


## Deflation
## when a Q[k] matrix become a "P" matrix, we deflate. This is checked by the sine term being basically 0.
function check_deflation{T}(state::RealDoubleShift{T}, tol = eps(T))
    for k in state.ctrs.stop_index:-1:state.ctrs.start_index
        if abs(vals(state.Q[k])[2]) <= tol
            deflate(state, k)
            return
        end
    end
end

# deflate a term
# turn on `show_status` to view sequence
function deflate{T}(state::RealDoubleShift{T}, k)

    # make a P matrix
    vals!(state.Q[k], getp(state.Q[k]), zero(T))
    
    # shift zero counter
    state.ctrs.zero_index = k      # points to a matrix Q[k] either RealRotator(-1, 0) or RealRotator(1, 0)
    state.ctrs.start_index = k + 1

    # reset counter
    state.ctrs.it_count = 1
end
