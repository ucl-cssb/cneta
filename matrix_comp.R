## Problem
# We want to find p(X -> Y)
# S1 and S2 are copy number profiles
# v1 v2 v3 v4 .. vn span the vector space and represent the possible mutational events
# we need to solve S2 - S1 = x1 v1 + x2 v2 + x3 v3 +  ... xn vn, where x1, x2.. xn are scalars
# can write this as:
#
#      [v1 v2 vn]x = y = S2 - S1
#
# so the problem is in general underdetermined
#
# [v1 v2 .. vn] in R(m x n). In general m < n since m is number of loci and n number of possible copy number changes
# x is underspecified: many choices of x lead to the same y


## Numerical example
# number of points
M = 10

# change basis
v1   = c(1,0,0,0,0,0,0,0,0,0)
v2   = c(0,1,0,0,0,0,0,0,0,0)
v3   = c(0,0,1,0,0,0,0,0,0,0)
v4   = c(0,0,0,1,0,0,0,0,0,0)
v5   = c(0,0,0,0,1,0,0,0,0,0)
v6   = c(0,0,0,0,0,1,0,0,0,0)
v7   = c(0,0,0,0,0,0,1,0,0,0)
v8   = c(0,0,0,0,0,0,0,1,0,0)
v9   = c(0,0,0,0,0,0,0,0,1,0)
v10  = c(0,0,0,0,0,0,0,0,0,1)
v11  = c(1,1,1,1,1,1,1,1,1,1)


S2 = c(0,0,0,0,0,1,0,1,0,1)
S1 = c(0,0,0,0,0,0,0,0,0,1)

V = as.matrix( cbind( v1, v2, v3, v4, v5, v6, v7, v8, v9, v10, v11 ) )

# the answer here is: y = 1*v6 + 1*v8
xsol = c(0,0,0,0,0,1,0,1,0,0,0)

# can check with
ysol = V %*% xsol

## Least-norm solution is given by
xln = t(V) %*% solve( V %*% t(V) ) %*% (S2-S1) 