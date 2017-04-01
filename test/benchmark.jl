using AMVW
using BenchmarkTools
using Polynomials
using PolynomialRoots
using DataFrames

function _residual_check(rs, rts)
    p = poly(rs)
    pp = polyder(p)

    # r1 = |P(lambda)/P'(lambda)|
    # r2 = |P(lambda)/P'(lambda)/lambda|

    r0 = maximum(norm.(sort(rs, by=norm) - sort(rts, by=norm)))  # nonsensical for some cases
    r1 = maximum(norm.(p.(rts) ./ pp.(rts)))
    r2 = maximum(norm.(p.(rts) ./ pp.(rts) ./ rts ))

    [r0, r1, r2]
end

function residual_check(rs)
    A = DataFrame(Polynomials=zeros(3),AMVW=zeros(3),PolynomialRoots=zeros(3))

    p = poly(rs)
    A[:,1] = _residual_check(rs, Polynomials.roots(p))
    A[:,2] = _residual_check(rs, AMVW.poly_roots(p.a))
    A[:,3] = _residual_check(rs, PolynomialRoots.roots(p.a))

    A
end


# small n, real
n = 10
rs = linspace(1/n, 1, n)
p = poly(rs)

@benchmark Polynomials.roots(p)
# julia> @benchmark Polynomials.roots(p)
# BenchmarkTools.Trial: 
#   memory estimate:  8.20 KiB
#   allocs estimate:  50
#   --------------
#   minimum time:     98.108 μs (0.00% GC)
#   median time:      102.220 μs (0.00% GC)
#   mean time:        108.525 μs (1.16% GC)
#   maximum time:     3.642 ms (94.16% GC)
#   --------------
#   samples:          10000
#   evals/sample:     1
#   time tolerance:   5.00%
#   memory tolerance: 1.00%

@benchmark AMVW.poly_roots(p.a)
# julia> @benchmark AMVW.poly_roots(p.a)
# BenchmarkTools.Trial: 
#   memory estimate:  4.27 KiB
#   allocs estimate:  66
#   --------------
#   minimum time:     93.309 μs (0.00% GC)
#   median time:      150.606 μs (0.00% GC)
#   mean time:        155.830 μs (0.26% GC)
#   maximum time:     4.237 ms (94.03% GC)
#   --------------
#   samples:          10000
#   evals/sample:     1
#   time tolerance:   5.00%
#   memory tolerance: 1.00%




@benchmark PolynomialRoots.roots(p.a)

# with deprecation warnings
# BenchmarkTools.Trial: 
#   memory estimate:  12.20 KiB
#   allocs estimate:  77
#   --------------
#   minimum time:     229.386 μs (0.00% GC)
#   median time:      241.773 μs (0.00% GC)
#   mean time:        253.786 μs (0.62% GC)
#   maximum time:     8.150 ms (0.00% GC)
#   --------------
#   samples:          10000
#   evals/sample:     1
#   time tolerance:   5.00%
#   memory tolerance: 1.00%

# Complex polynomials
rs = [x+im for x in 1.0 : 6]
p = poly(rs)

@benchmark Polynomials.roots(p)

# julia> @benchmark Polynomials.roots(p)
# BenchmarkTools.Trial: 
#   memory estimate:  5.59 KiB
#   allocs estimate:  52
#   --------------
#   minimum time:     33.073 μs (0.00% GC)
#   median time:      34.435 μs (0.00% GC)
#   mean time:        36.119 μs (1.96% GC)
#   maximum time:     3.713 ms (95.31% GC)
#   --------------
#   samples:          10000
#   evals/sample:     1
#   time tolerance:   5.00%
#   memory tolerance: 1.00%

@benchmark AMVW.poly_roots(p.a)
# BenchmarkTools.Trial: 
#   memory estimate:  4.66 KiB
#   allocs estimate:  56
#   --------------
#   minimum time:     102.511 μs (0.00% GC)
#   median time:      143.153 μs (0.00% GC)
#   mean time:        157.052 μs (0.70% GC)
#   maximum time:     49.993 ms (0.00% GC)
#   --------------
#   samples:          10000
#   evals/sample:     1
#   time tolerance:   5.00%
# memory tolerance: 1.00%

@benchmark PolynomialRoots.roots(p.a)
# BenchmarkTools.Trial: 
#   memory estimate:  11.63 KiB
#   allocs estimate:  76
#   --------------
#   minimum time:     226.424 μs (0.00% GC)
#   median time:      232.025 μs (0.00% GC)
#   mean time:        248.389 μs (0.89% GC)
#   maximum time:     7.717 ms (0.00% GC)
#   --------------
#   samples:          10000
#   evals/sample:     1
#   time tolerance:   5.00%
# memory tolerance: 1.00%

## XXX We have an issue with this one though
n = 10
ts = linspace(1/n, 1.0, n) * 2pi
rs = [complex(cos(t), sin(t)) for t in ts]
p = poly(rs) 

@benchmark Polynomials.roots(p)
# julia> @benchmark PolynomialRoots.roots(p.a)
# BenchmarkTools.Trial: 
#   memory estimate:  11.95 KiB
#   allocs estimate:  76
#   --------------
#   minimum time:     301.440 μs (0.00% GC)
#   median time:      306.313 μs (0.00% GC)
#   mean time:        326.734 μs (0.61% GC)
#   maximum time:     4.879 ms (91.35% GC)
#   --------------
#   samples:          10000
#   evals/sample:     1
#   time tolerance:   5.00%
#   memory tolerance: 1.00%

@benchmark AMVW.poly_roots(p.a) #  use state.ray=false
# ## XXX Stil misses alot
# BenchmarkTools.Trial: 
#   memory estimate:  5.88 KiB
#   allocs estimate:  68
#   --------------
#   minimum time:     354.586 μs (0.00% GC)
#   median time:      401.278 μs (0.00% GC)
#   mean time:        443.340 μs (0.30% GC)
#   maximum time:     6.563 ms (90.24% GC)
#   --------------
#   samples:          10000
#   evals/sample:     1
#   time tolerance:   5.00%
#   memory tolerance: 1.00%


@benchmark PolynomialRoots.roots(p.a)

# julia> @benchmark PolynomialRoots.roots(p.a)
# BenchmarkTools.Trial: 
#   memory estimate:  11.95 KiB
#   allocs estimate:  76
#   --------------
#   minimum time:     302.533 μs (0.00% GC)
#   median time:      309.085 μs (0.00% GC)
#   mean time:        330.070 μs (0.60% GC)
#   maximum time:     5.024 ms (88.43% GC)
#   --------------
#   samples:          10000
#   evals/sample:     1
#   time tolerance:   5.00%
#   memory tolerance: 1.00%

# "Big" polynomials
n = 10
rs = linspace(1/n, big(1.0), n)

p = poly(rs) 

@benchmark Polynomials.roots(p)
## Error (no support for big)

@benchmark AMVW.poly_roots(p.a)
# julia> @benchmark AMVW.poly_roots(p.a)
# BenchmarkTools.Trial: 
#   memory estimate:  10.14 MiB
#   allocs estimate:  204777
#   --------------
#   minimum time:     26.071 ms (0.00% GC)
#   median time:      56.661 ms (33.13% GC)
#   mean time:        52.697 ms (24.94% GC)
#   maximum time:     74.463 ms (25.16% GC)
#   --------------
#   samples:          95
#   evals/sample:     1
#   time tolerance:   5.00%
#   memory tolerance: 1.00%

# julia> maximum(norm.(sort(AMVW.poly_roots(p.a), by=norm) - sort(rs, by=norm)))
# 3.542606431077360733112308774961528146389295822293876354568756300943283335489686e-71

@benchmark PolynomialRoots.roots(p.a) 
# BenchmarkTools.Trial: 
#   memory estimate:  1.25 MiB
#   allocs estimate:  26226
#   --------------
#   minimum time:     1.669 ms (0.00% GC)
#   median time:      1.873 ms (0.00% GC)
#   mean time:        3.034 ms (36.15% GC)
#   maximum time:     26.625 ms (91.54% GC)
#   --------------
#   samples:          1621
#   evals/sample:     1
#   time tolerance:   5.00%
#   memory tolerance: 1.00%

# julia> maximum(norm.(sort(PolynomialRoots.roots(p.a), by=norm) - sort(rs, by=norm)))
# 1.154958001715555551976043766371611436643198308737886091026144983582996624795927e-72

##################################################

# larger n (50)

n = 50
rs = linspace(1/n, 1, n)
p = poly(rs)

@benchmark Polynomials.roots(p)
# julia> @benchmark Polynomials.roots(p)
# BenchmarkTools.Trial: 
#   memory estimate:  61.78 KiB
#   allocs estimate:  82
#   --------------
#   minimum time:     1.006 ms (0.00% GC)
#   median time:      1.073 ms (0.00% GC)
#   mean time:        1.128 ms (0.57% GC)
#   maximum time:     4.822 ms (65.08% GC)
#   --------------
#   samples:          4306
#   evals/sample:     1
#   time tolerance:   5.00%
#   memory tolerance: 1.00%


@benchmark AMVW.poly_roots(p.a)
# julia> @benchmark AMVW.poly_roots(p.a)
# BenchmarkTools.Trial: 
#   memory estimate:  12.73 KiB
#   allocs estimate:  186
#   --------------
#   minimum time:     1.556 ms (0.00% GC)
#   median time:      2.432 ms (0.00% GC)
#   mean time:        2.515 ms (0.08% GC)
#   maximum time:     8.164 ms (46.20% GC)
#   --------------
#   samples:          1962
#   evals/sample:     1
#   time tolerance:   5.00%
#   memory tolerance: 1.00%

@benchmark PolynomialRoots.roots(p.a)
# julia> @benchmark PolynomialRoots.roots(p.a)
# BenchmarkTools.Trial: 
#   memory estimate:  14.72 KiB
#   allocs estimate:  77
#   --------------
#   minimum time:     422.714 μs (0.00% GC)
#   median time:      426.969 μs (0.00% GC)
#   mean time:        446.475 μs (0.53% GC)
#   maximum time:     4.512 ms (88.61% GC)
#   --------------
#   samples:          10000
#   evals/sample:     1
#   time tolerance:   5.00%
#   memory tolerance: 1.00%

##################################################

# residual check
n = 10;
rs = linspace(1/n, 1, n);
p = poly(rs);
residual_check(rs)

# julia> residual_check(rs)
# 3×3 DataFrames.DataFrame
# │ Row │ Polynomials │ AMVW       │ PolynomialRoots │
# ├─────┼─────────────┼────────────┼─────────────────┤
# │ 1   │ 2.09855e-10 │ 1.18857e-9 │ 5.75882e-11     │
# │ 2   │ 1.4072e-10  │ 1.13841e-9 │ 3.38759e-11     │
# │ 3   │ 2.01029e-10 │ 2.06593e-9 │ 4.23449e-11     │

n = 15;
rs = linspace(1/n, 1, n);
p = poly(rs);
residual_check(rs)

# julia> residual_check(rs)
# 3×3 DataFrames.DataFrame
# │ Row │ Polynomials │ AMVW       │ PolynomialRoots │
# ├─────┼─────────────┼────────────┼─────────────────┤
# │ 1   │ 2.23049e-6  │ 2.742e-5   │ 1.25911e-7      │
# │ 2   │ 2.30464e-6  │ 2.74255e-5 │ 2.07169e-7      │
# │ 3   │ 3.14269e-6  │ 4.83827e-5 │ 3.10753e-7      │


n = 50;
rs = linspace(1/n, 1, n);
p = poly(rs);
residual_check(rs)

# julia> residual_check(rs)
# 3×3 DataFrames.DataFrame
# │ Row │ Polynomials │ AMVW       │ PolynomialRoots │
# ├─────┼─────────────┼────────────┼─────────────────┤
# │ 1   │ 0.943682    │ 0.640697   │ 0.594116        │
# │ 2   │ 0.0213801   │ 0.0183637  │ 0.0446998       │
# │ 3   │ NaN         │ 1.93315e10 │ 0.0488876       │

n = 10
ts = linspace(1/n, 1.0, n) * 2pi
rs = [complex(cos(t), sin(t)) for t in ts]
p = poly(rs) 
residual_check(rs)

# 3×3 DataFrames.DataFrame
# │ Row │ Polynomials │ AMVW        │ PolynomialRoots │
# ├─────┼─────────────┼─────────────┼─────────────────┤
# │ 1   │ 1.90211     │ 1.90211     │ 1.90211         │
# │ 2   │ 3.39034e-15 │ 3.52069e-15 │ 3.5958e-16      │
# │ 3   │ 3.39034e-15 │ 3.52069e-15 │ 3.5958e-16      │

