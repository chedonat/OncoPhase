library(limSolve)

# predefined some scenario and observed VAF
M = matrix( c(1, 1, 2, 0, 1, 2), nrow=2, byrow=TRUE )
C = matrix( c(2, 2, 3, 2, 2, 3), nrow=2, byrow=TRUE )
w = matrix( c(0.584, 0, 0, 0.294), nrow=2, byrow=TRUE )

# main linear system
A = w%*%C-M
b = matrix( c(0, 0), nrow=2, byrow=TRUE)

# equality constraints
e = matrix( c(1, 1, 1), nrow=1, byrow=TRUE) # theta_G + theta_A + theta_B = 1
f = 1

# bounds (first 3 rows theta_G, theta_A, theta_B > 0, second 3 rows -theta_G, -theta_A, -theta_B > -1 (equiv. theta_G, theta_A, theta_b < 1))
g = matrix(c(1, 0, 0, 0, 1, 0, 0, 0, 1, -1, 0, 0, 0, -1, 0, 0, 0, -1), nrow=6, byrow=TRUE)
h = matrix(c(0, 0, 0, -1, -1, -1), nrow=6, byrow=TRUE)

# use constrained linear solver
output_without_bounds = lsei( A = A, B = b, E = e, F = f)
output_with_bounds = lsei( A = A, B = b, E = e, F = f, G = g, H = h)

print(output_without_bounds$X)
print(output_with_bounds$X)




