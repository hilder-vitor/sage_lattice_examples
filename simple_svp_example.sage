from sage.modules.free_module_integer import IntegerLattice as Lattice

n = 10

B = 5*Matrix.random(ZZ, n, n) + 2*Matrix.random(ZZ, n, n)

print "Basis of L:"
print B

# create a lattice with basis given by the rows of B
L = Lattice(B)

det = L.volume()

v = L.shortest_vector()

l1 = v.norm()
print "lambda_1(L) =", l1.n()
print "Shortest vector:", v
# assert Minkowski theorem is valid
assert (l1 <= sqrt(n) * (det ** (1/n)))
print "lambda_1(L) <= sqrt(n) * (det(L) ^ (1/n)) ? True"

