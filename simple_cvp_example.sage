from sage.modules.free_module_integer import IntegerLattice as Lattice

n = 4

B = Matrix.random(ZZ, n, n)

print "Basis of L:"
print B

# create a lattice with basis given by the rows of B
L = Lattice(B, lll_reduce=True)

target = vector(ZZ, [0]*n)

# selecting a target not in L, otherwise it will be the closest vector itself
while target in L:
    u = L.random_element()
    target = u + vector(ZZ, [0]*(n-1) + [1])

print "u =", u
print "target =", target

v = L.closest_vector(target)

print "CVP(target) =", v
print "CVP(target) == u? ", (v == u)

