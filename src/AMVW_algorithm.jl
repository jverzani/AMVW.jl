
## Main algorithm of AMV&W
## This follows that given in the paper very closely
function AMVW_algorithm{T,St,P,Tw}(state::FactorizationType{T, St, P, Tw})


    it_max = 20 * state.N
    kk = 0

    while kk <= it_max

        ## finished up!
        state.ctrs.stop_index <= 0 && return

        state.ctrs.it_count += 1

        
        check_deflation(state)
        kk += 1

#        show_status(state)

        k = state.ctrs.stop_index

        if state.ctrs.stop_index - state.ctrs.zero_index >= 2
            
            bulge_step(state)
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



## We decompose ps into qs or vs, ws for pencil
## caller might look like
function basic_decompose{T}(ps::Vector{T})
#    ps, k = deflate_leading_zeros(ps)  ## ASSUMED
    qs = ps[1:end-1] / ps[end]
    par = iseven(length(qs)) ? one(T) : -one(T)
    qs = vcat(-qs[2:end], par*qs[1], -par * one(T))
    qs
end

function basic_decompose1(ps)
#    ps, k = deflate_leading_zeros(ps)  ## ASSUMED
    qs = reverse_poly(ps)  #
    qs
end


function basic_decompose_1(ps)
    #    ps, k = deflate_leading_zeros(ps)  ## ASSUMED
    qs = ps[1:end-1] / ps[end]
    head = shift!(qs)
    push!(qs, head)
    qs
end

function basic_decompose_pencil{T}(ps::Vector{T})
    N = length(ps)-1
    qs = zeros(T, N+1)
    qs[N]  = ps[end]
    qs[N+1] = -one(T)

    par = iseven(N) ? one(T) : -one(T)
    ps = vcat(-ps[2:end-1], par * ps[1], -par*one(T))

    ps, qs
    
end

## for real like 8e-6 * N^2 run time
##     allocations like c + 3n where c covers others.
## for complex   1e-5 * N^2 run time (1.25 more)
##     allocations like  3n too?
function amvw{T <: Real}(ps::Vector{T})
    state = convert(FactorizationType{T, Val{:DoubleShift}, Val{:NoPencil}, Val{:NotTwisted}}, ps)
    init_state(state, basic_decompose)
    state
end

function amvw{T <: Real}(ps::Vector{Complex{T}})
    state = convert(FactorizationType{T, Val{:SingleShift}, Val{:NoPencil}, Val{:NotTwisted}}, ps)
    init_state(state, basic_decompose)
    state
end

function amvw_pencil{T <: Real}(ps::Vector{T}, decompose=basic_decompose_pencil)
    state = convert(FactorizationType{T, Val{:DoubleShift}, Val{:HasPencil}, Val{:NotTwisted}}, ps)
    init_state(state, decompose)
    state
end

function amvw_pencil{T <: Real}(ps::Vector{Complex{T}}, decompose=basic_decompose_pencil)
    state = convert(FactorizationType{T, Val{:SingleShift}, Val{:HasPencil}, Val{:NotTwisted}}, ps)
    init_state(state, decompose)
    state
end



function poly_roots{T}(ps::Vector{T})
    # deflate 0s
    ps, K = deflate_leading_zeros(ps)

    state = amvw(ps)
    AMVW_algorithm(state)

    if K > 0
        ZS = zeros(T, K)
        if T <: Complex
            append!(state.REIGS, real.(ZS))
            append!(state.IEIGS, real.(ZS))
        else
            append!(state.REIGS, ZS)
            append!(state.IEIGS, ZS)
        end            
    end

    complex.(state.REIGS, state.IEIGS)
end
