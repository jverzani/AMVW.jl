#__precompile__(true)
module AMVW

using Compat
# package code goes here


include("types.jl")
include("utils.jl")
include("transformations.jl")
include("factorization.jl")
include("diagonal-block.jl")
include("bulge.jl")
include("deflation.jl")
include("AMVW_algorithm.jl")

include("diagnostics.jl")

end # module
