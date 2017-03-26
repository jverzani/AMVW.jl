## Types
@compat abstract type CoreTransform{T} end
@compat abstract type Rotator{T} <: CoreTransform{T} end


#the index is supeflous for now, and a bit of a hassle to keep immutable
#but might be of help later if twisting is approached. Shouldn't effect speed, but does mean 9N storage, not 6N
#so may be
#
mutable struct RealRotator{T} <: Rotator{T}
c::T
s::T
i::Int
end

function Base.ctranspose(r::RealRotator)
    RealRotator(r.c, -r.s, r.i)
end


Base.one{T}(::Type{RealRotator{T}}) = RealRotator(one(T), zero(T), 0)
Base.ones{T}(S::Type{RealRotator{T}}, N) = [one(S) for i in 1:N]

## get/set values
vals{T}(r::RealRotator{T}) = (r.c,r.s)
function vals!{T}(r::RealRotator, c::T, s::T)
    # normalize in case of roundoff errors
    # but, using hueristic on 6.3 on square roots
    
    nrmi = c^2 + s^2 
    nrmi = norm(nrmi - one(T)) >= 1e2*eps(T) ? inv(sqrt(nrmi)) : one(T)
    r.c = c * nrmi
    r.s = s * nrmi
end
idx(r::RealRotator) = r.i
idx!(r::RealRotator, i::Int) = r.i = i

Base.copy(a::RealRotator) = RealRotator(a.c, a.s, a.i) #a.c, a.s, a.i)
function Base.copy!(a::RealRotator, b::RealRotator)
    vals!(a, vals(b)...)
    idx!(a, idx(b))
end
    

# Core transform is 2x2 matrix [a b; c d]
mutable struct  RealTransform{T} <: CoreTransform{T}
    xs::Vector{T} # [a b; c d]
    i::Int
end
Base.ctranspose(r::RealTransform) = RealTransform(r.xs[[1,3,2,4]], r.i)


## A container for our counters
mutable struct DoubleShiftCounter
    zero_index::Int
    start_index::Int
    stop_index::Int
    it_count::Int
    tr::Int
end

@compat abstract type ShiftType{T} end
struct RealDoubleShift{T} <: ShiftType{T} 
    N::Int
    POLY::Vector{T}
    Q::Vector{RealRotator{T}}
    Ct::Vector{RealRotator{T}}  # We use C', not C here
    B::Vector{RealRotator{T}}
    REIGS::Vector{T}
    IEIGS::Vector{T}
    ## reusable storage
U::RealRotator{T}
Ut::RealRotator{T}
V::RealRotator{T}
Vt::RealRotator{T}
    W::RealRotator{T}
    A::Matrix{T}    # for parts of A = QR
    Bk::Matrix{T}   # for diagonal block
    R::Matrix{T}    # temp storage, sometimes R part of QR
    e1::Vector{T}   # eigen values e1, e2
    e2::Vector{T}
    ctrs::DoubleShiftCounter
end

function Base.convert{T}(::Type{RealDoubleShift}, ps::Vector{T})
    N = length(ps)
    
    RealDoubleShift(N, ps,
                    ones(RealRotator{T}, N), #Q
                    ones(RealRotator{T}, N), #Ct
                    ones(RealRotator{T}, N), #B
                    zeros(T, N),  zeros(T, N), #EIGS
                    one(RealRotator{T}), one(RealRotator{T}),
                    one(RealRotator{T}), one(RealRotator{T}),
                    one(RealRotator{T}), #U,U',V,V',W
                    zeros(T, 2, 2),zeros(T, 3, 2),zeros(T, 3, 2), # A Bk R
                    zeros(T,2), zeros(T,2),
                    DoubleShiftCounter(0,1,N-1, 0, N-2)
    )
end

##################################################
## We use two complex, rather than 3 reals here.
## Same storage, as we don't need to include a D
    
mutable struct ComplexRotator{T} <: Rotator{T}
c::Complex{T}
s::Complex{T}
i::Int
end

function Base.ctranspose(r::ComplexRotator)
    ComplexRotator(conj(r.c), -r.s, r.i)
end


Base.one{T}(::Type{ComplexRotator{T}}) = ComplexRotator(complex(one(T), zero(T)), complex(zero(T), zero(T)), 0)
Base.ones{T}(S::Type{ComplexRotator{T}}, N) = [one(S) for i in 1:N]

## get/set values
vals{T}(r::ComplexRotator{T}) = (r.c, r.s)
function vals!{T}(r::ComplexRotator, c::Complex{T}, s::Complex{T})
    # normalize in case of roundoff errors
    # but, using hueristic on 6.3 on square roots
    
    nrmi = sqrt(abs(c * conj(c) + s * conj(s)))
    nrmi = norm(nrmi - one(T)) >= eps(T) ? inv(sqrt(nrmi)) : one(T)
    r.c = c * nrmi
    r.s = s * nrmi
end
vals!{T}(r::ComplexRotator, c::Complex{T}, s::T) = vals!(r, c, complex(s,zero(T)))
vals!{T}(r::ComplexRotator{T}, c::T, s::T) = vals!(r, complex(c, zero(T)), complex(s, zero(T)))
idx(r::ComplexRotator) = r.i
idx!(r::ComplexRotator, i::Int) = r.i = i

Base.copy(a::ComplexRotator) = ComplexRotator(a.c, a.s, a.i)
function Base.copy!(a::ComplexRotator, b::ComplexRotator)
    vals!(a, vals(b)...)
    idx!(a, idx(b))
end
  

## A container for our counters
## XXX This is not complex/real dependent
## XXX part of algorithm
mutable struct SingleShiftCounter
    zero_index::Int
    start_index::Int
    stop_index::Int
    it_count::Int
    tr::Int
end

@compat abstract type ShiftType{T} end
struct ComplexSingleShift{T} <: ShiftType{T} 
    N::Int
    POLY::Vector{Complex{T}}
    Q::Vector{ComplexRotator{T}}
    Ct::Vector{ComplexRotator{T}}  # We use C', not C here
    B::Vector{ComplexRotator{T}}  
    REIGS::Vector{T}
    IEIGS::Vector{T}
    ## reusable storage
U::ComplexRotator{T}
Ut::ComplexRotator{T}
    A::Matrix{Complex{T}}    # for parts of A = QR
    Bk::Matrix{Complex{T}}   # for diagonal block
    R::Matrix{Complex{T}}    # temp storage, sometimes R part of QR
    e1::Vector{T}   # eigen values e1, e2, store as (re,imag)
e2::Vector{T}
ray::Bool
    ctrs::SingleShiftCounter
end

function Base.convert{T}(::Type{ComplexSingleShift}, ps::Vector{Complex{T}})
    N = length(ps)
    
    ComplexSingleShift(N, ps,
                       ones(ComplexRotator{T}, N), #Q
                       ones(ComplexRotator{T}, N), #Ct
                       ones(ComplexRotator{T}, N), #B
                       zeros(T, N),  zeros(T, N), #EIGS
                       one(ComplexRotator{T}), one(ComplexRotator{T}), #U, Ut
                       zeros(Complex{T}, 2, 2),zeros(Complex{T}, 3, 2),
                       zeros(Complex{T}, 3, 2), # A Bk R
    zeros(T,2), zeros(T,2),
    true,  # true for Wilkinson, 1 for Rayleigh
    SingleShiftCounter(0,1,N-1, 0, N-2)
    )
end

