#__precompile__(true)
module AMVW

using Compat
# package code goes here


include("types.jl")
include("utils.jl")
include("transformations.jl")
include("bulge.jl")
include("factorization.jl")
include("diagnostics.jl")
include("AMVW_algorithm.jl")


end # module
