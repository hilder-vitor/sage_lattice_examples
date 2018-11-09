from sage.modules.free_module_integer import IntegerLattice as Lattice

n = 7

# Create the knapsack instance
density = 3/n  # must be in the interval ]0, 1]
upper_bound_a = 2**(int(ceil(n / density)))
while True:
    a = vector([randint(1, upper_bound_a)  for i in xrange(n)])
    x = vector([randint(0, 1) for i in xrange(n)]) # solution
    s = a * x
    if a != 0 and x != 0 and s != 0:
        break
print "a =", a
print "x =", x
print "s =", s

# construct the basis with -s and a in the first column.
# For instance, for n=3, B will be as follows
#   [ -s 0 0 0 ]
#   [ a0 1 0 0 ]
#   [ a1 0 1 0 ]
#   [ a2 0 0 1 ]
B = block_matrix(ZZ, 2, 2, [[-s, 0], [Matrix(ZZ, n, 1, list(a)), 1]])

print "Basis of L:"
print B

# create a lattice with basis given by the rows of B
L = Lattice(B)

L.LLL()
v = L.reduced_basis[0] # select the first vector of the reduced basis as a solution
solution = v[1:]

print "Found solution:", solution
print "Is x == found solution?", x == solution
print "is (found solution)*a = s?", (solution * a == s)
