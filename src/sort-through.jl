module T2

using AMVW


## XXX This is just junk for now to sort through

## How to set up shifttype
# ShiftType{T, ShiftType( :SS, :DS), PencilType ("Pencil, NoPencil), Twisted ("Twisted", NotTwisted)}(p, args...)
abstract type ShiftType{T, ST} end

struct Real_DoubleShift_NoPencil_NotTwisted{T} <: ShiftType{T, Val{:DS}}
U::T
V::T
Q::Vector{T}
C::Vector{T}
B::Vector{T}
end


struct Complex_SingleShift_NoPencil_NotTwisted{T} <: ShiftType{T, Val{:SS}}
U::T
V::T
D::Vector{T}
Q::Vector{T}
C::Vector{T}
B::Vector{T}
end

function Base.convert{T}(::Type{ShiftType{T,Val{:DS}}}, a::Vector{T})
    RealDoubleShift(1.0, 1.0, a,a,a)
end

function Base.convert{T}(::Type{ShiftType{T,Val{:SS}}}, a::Vector{T})
    ComplexSingleShift(1.0, 1.0,a, a,a,a)
end

f{T}(a::ShiftType{T, Val{:DS}}) = "DS"
f{T}(a::ShiftType{T, Val{:SS}}) = "SS"

g{T,S}(a::ShiftType{T, S}) = "match"





vals = AMVW.vals
vals! = AMVW.vals!
idx = AMVW.idx
idx! = AMVW.idx!
givensrot = AMVW.givensrot
reverse_poly = AMVW.reverse_poly
as_full = AMVW.as_full
ShiftType = ComplexRealSingleShift = AMVW.ComplexRealSingleShift
CRR = ComplexRealRotator = AMVW.ComplexRealRotator

# pass through a upper triangle
# matrix described by sequenc so Cs going up, Bs down
# T U --> U T
function through_triu(U::Rotator, Cs, Bs, ::Type{Val{:right}})
    i = idx(U)
    turnover(Bs[i], Bs[i+1], U, Val{:right})
    turnover(Cs[i+1], Cs[i], U, Val{:right})
end
# U T --> T U
function through_triu(U::Rotator, Cs, Bs, ::Type{Val{:left}})
    i = idx(U)
    turnover(U, Cs[i+1], Cs[i], Val{:left})
    turnover(U, Bs[i], Bs[i+1], Val{:left})
end
# T^(-1) U --> U T^(-1)
function through_triu(U::Rotator, Cs, Bs, ::Type{Val{:right}}, ::Type{Val{:inv}})
    i = idx(U)
    copy!(U, conj(U))
    copy!(Bs[i], conj(Bs[i]))
    copy!(Bs[i+1], conj(Bs[i+1]))
    copy!(Cs[i], conj(Cs[i]))
    copy!(Cs[i+1], conj(Cs[i+1]))    
    # 
    turnover(U, Bs[i+1], Bs[i], Val{:left})
    turnover(U, Cs[i], Cs[i+1], Val{:left})
    copy!(U, conj(U))
    copy!(Bs[i], conj(Bs[i]))
    copy!(Bs[i+1], conj(Bs[i+1]))
    copy!(Cs[i], conj(Cs[i]))
    copy!(Cs[i+1], conj(Cs[i+1]))  
end

# U T^(-1) --> T^(-1) U
function through_triu(U::Rotator, Cs, Bs, ::Type{Val{:left}}, ::Type{Val{:inv}})
    i = idx(U)
    copy!(U, conj(U))
    copy!(Bs[i], conj(Bs[i]))
    copy!(Bs[i+1], conj(Bs[i+1]))
    copy!(Cs[i], conj(Cs[i]))
    copy!(Cs[i+1], conj(Cs[i+1]))    
    # 
    turnover(Cs[i], Cs[i+1], U, Val{:right})
    turnover(Bs[i+1], Bs[i], U, Val{:right})
    copy!(U, conj(U))
    copy!(Bs[i], conj(Bs[i]))
    copy!(Bs[i+1], conj(Bs[i+1]))
    copy!(Cs[i], conj(Cs[i]))
    copy!(Cs[i+1], conj(Cs[i+1]))    
end
    
# absorb Ut Q
# might need to account for phase (D)
# might need to account for twisting (sigma)
function through_Q(Ut, Q, ::Type{Val{:left}})
    i = idx(Ut)
    i > 1 && dflip(Ut, Q[i-1])
    fuse(Ut, Q[i], Val{:right})
end
function through_Q(Ut, Q, D::Vector{Complex{T}}, ::Type{Val{:left}})
    i = idx(Ut)
    alpha = fuse(Ut, Q[i], Val{:left})
    # absorb alpha, here we cascade but don't have stop index as a guide
    j = i + 1
    while j <= length(Q) && !is_identity(Q[j])
        c, s = vals(Q[j])
        vals!(Q[j], c*conj(alpha), s)
        j = j + 1
    end
    D[i] = alpha; D[j] = conj(alpha)
end

function through_Q(Ut, Q, sigma::Vector{Int}, ::Type{Val{:left}})
    ## This is tricky one!
end

function through_Q(Ut, Q, D, sigma, ::Type{Val{:left}})
    ## This is tricky one!
end

## from right
## here we return true if absorbed
## false if not
function through_Q(U, Q, ::Type{Val{:right}})
    i = idx(Ut)
    if i < length(Q) && !is_identity(Q[i+1])
        turnover(Q[i], Q[i+1], U, Val{:right})
        false
    else
        fuse(Q[i], U, Val{:left})
        true
    end
end

function through_Q(U, Q, D::Vector{Complex{T}}, ::Type{Val{:right}})
    i = idx(Ut)
    passthrough(view(D,i:i+1), U) # allocate

    if i < length(Q) && !is_identity(Q[i+1])
        turnover(Q[i], Q[i+1], U, Val{:right})
        false
    else
        alpha = fuse(Q[i], U, Val{:left})
        D[i] *= alpha
        D[i+1] *= conj(alpha)
        true
    end
end

function through_Q(Ut, Q, sigma::Vector{Int}, ::Type{Val{:right}})
    ## Twisted. This is tricky one!
end

function through_Q(Ut, Q, D, sigma, ::Type{Val{:right}})
    ## Twisted, This is tricky one!
end

# basic steps
# create_bulge
# through_Q(state.U, state.Q, [D], [sigma], Val{:left})
# flag = false
# while !flag
#   [through_triu(state.U, state.Ct1, state.B1, Val{:right})]
#   through_triu(state.U, state.Ct, state.B, Val{:right})
#   flag = through_Q(state.U, state.Q, [D], [sigma], Val{:right})
# end


## We have three things:
# * real double shift / complex single shift (RDS -- no D; CSS -- has D)
# * pencil or not (no pencil only Ct, B; has pencil Ct, B and Ct1, B1)
# * twisted, not twisted (not twisted, sigma[i] = i; twisted we have a sigma) Could also have cdc -- alternating)

# our differences:
# create bulge (depends on RDS or CSS)
# diagonal block (depends on pencil, RDS/CSS, twisted, not twisted) -- need to break up
# absorb Ut (Ut, Vt) depends on RDS, CSS; twisted or not;
# absort U (U,V) depends on RDS/CSS; pencil or not; depends on through Q

# three main things
## Ut (Q D CB) -> (Q D C B)
# absorb_transpose(state, RDS/CSS, Twisted/Not)
function absorb_transpose(state, ::Type{Val{:RDS}}, ::Type{Val{:NotTwisted}})
    # absorb Ut, Vt, leave W
    # not twisted
end

function absorb_transpose(state, ::Type{Val{:RDS}}, ::Type{Val{:Twisted}})
    # absorb Ut, Vt, leave W
    # not twisted
end

function absorb_transpose(state, ::Type{Val{:CSS}}, ::Type{Val{:NotTwisted}})
    through_Q(state.Ut, state.Q, state.D, Val{:left})
end

function absorb_transpose(state, ::Type{Val{:CSS}}, ::Type{Val{:Twisted}})
    through_Q(state.Ut, state.Q, state.D, state.sigma, Val{:left})
end

## pass through triu part R or VW^{-1}
## need direction
function pass_triu(state, dir,  ::Type{Val{:RDS}}, ::Type{Val{:NonPencil}})
    ## end!
    through_triu(state.V, state.Ct, state.B, Val{:right})
    through_triu(state.U, state.Ct, state.B, Val{:right})        
end

function pass_triu(state,  ::Type{Val{:RDS}}, ::Type{Val{:Pencil}})
end

function pass_triu(state, dir,  ::Type{Val{:CSS}}, ::Type{Val{:NonPencil}})
    through_triu(state.U, state.Ct, state.B, dir)
end

function pass_triu(state,  dir, ::Type{Val{:CSS}}, ::Type{Val{:Pencil}})
    if dir == Val{:right}
        through_triu(state.U, state.Ct1, state.B1,dir, Val{:inv})    
        through_triu(state.U, state.Ct, state.B, dir)
    else
        through_triu(state.U, state.Ct, state.B, dir)
        through_triu(state.U, state.Ct1, state.B1,dir, Val{:inv})
    end
end




## bulge step for
## ComplexRealReal/NotPencil/NotTwisted
function one_bulge_step(state::CSS)
    AMVW.create_bulge(state)
    copy(state.Ut, conj(state.U))
    through_Q(state.Ut, state.Q, Val{:left})
    flag = false
    while !flag
        through_triu(state.U, state.Ct, state.B, Val{:right})
        flag = through_Q(state.U, state.Q, Val{:right})
    end
end
         
function one_bulge_step(state::RDS)
    AMVW.create_bulge(state)

    copy(state.Ut, conj(state.U))
    copy(state.Vt, conj(state.V))    

    # through Q leaves W
    i = idx(state.Ut)
    turnover(state.Ut, state.Vt, state.Q[i], Val{:right})
    copy!(state.W, state.Q[i])
    copy!(state.Q[i], state.Ut)
    fuse(state.Vt, state.Q[i+1], Val{:right})
    
    uflag = vflag = false
    while !uflag
        through_triu(state.V, state.Ct, state.B, Val{:right})
        through_triu(state.U, state.Ct, state.B, Val{:right})
        if !vflag
            vflag = through_Q(state.V, state.Q, Val{:right})
        end
        uflag = through_Q(state.U, state.Q, Val{:right})        
    end
end
   

function one_bulge_step(state::CSS, ::Type{Val{:Pencil}})
    AMVW.create_bulge(state)
    copy(state.Ut, conj(state.U))
    through_Q(state.Ut, state.Q, Val{:left})
    flag = false
    while !flag
        through_triu(state.U, state.Ct1, state.B1, Val{:right}, Val{:inv})
        through_triu(state.U, state.Ct, state.B, Val{:right})        
        flag = through_Q(state.U, state.Q, Val{:right})
    end
end
         
function one_bulge_step(state::RDS, ::Type{Val{:Pencil}})
    AMVW.create_bulge(state)
    copy!(state.Ut, conj(state.U))
    copy!(state.Vt, conj(state.V))

    # through Q leaves W
    i = idx(state.Ut)
    turnover(state.Ut, state.Vt, state.Q[i], Val{:right})
    copy!(state.W, state.Q[i])
    copy!(state.Q[i], state.Ut)
    fuse(state.Vt, state.Q[i+1], Val{:right})

    
    flag = false
    while !flag
        through_triu(state.V, state.Ct1, state.B1, Val{:right}, Val{:inv})
        through_triu(state.V, state.Ct, state.B, Val{:right})
        through_triu(state.U, state.Ct1, state.B1, Val{:right}, Val{:inv})
        through_triu(state.U, state.Ct, state.B, Val{:right})
        flag = through_Q(state.U, state.Q, Val{:right})
    end
end
   





# absorb (Q D) U or move on
# return true if absorbed, false if not


type TwistedQReal
    Q
    sigma
end
  
type TwistedQComplexReal
    Q
    sigma
    D
end
    



# bulge chasing
# A.create_bulge(state)

## left_Q
function left_Q(Ut, Qtype)
    # u -> Q1
    #       Q2
    #        Q3

    
    

