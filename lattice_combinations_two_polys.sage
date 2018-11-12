from sage.modules.free_module_integer import IntegerLattice as Lattice

N = 6  # Degree of polynomial used as modulus and dimension of lattice

q = random_prime(1000, 500)

R.<x> = ZZ['x']
f = x^N + 1     # modulus defining R / <f>

print "Ring: ZZ[x]/<", f, ">"
print "q =", q

def phi(h): # isomorphism from R/f to Z^n (return a vector representing h)
    v = h.coefficients(sparse=False)
    return v + [0]*(N - len(v))

def inv_phi(vec_h): # isomorphism from Z^n to R/f (return the polynomial represented by vec_h)
    return R(vec_h)


# Return a matrix H such that for any p in R/<f>, H*phi(p) equals h*p % f,
# that is, the matrix that represents the product by the polynomial h.
def matrix_polynomial_product(h):
    H = Matrix(ZZ, [phi(h * x^i % f) for i in xrange(N)])
    return H.transpose()

# All the vectors in L are supposed to have the format 
#     phi(u*q + c1*a % f + c2*b % f) || phi(a)|| phi(b)
# for some vectors u, a, and b of R/<f>.
# This function receives a v from L and return such polynomials u, a, and b.
def lattice_point_to_polys(v):
    v = list(v) # treat vector over ZZ as a simple list

    v0 = inv_phi(v[0:N]) # polynomial corresponding to the first N entries
    a = inv_phi(v[N:2*N])
    b = inv_phi(v[2*N:3*N])
    u = (v0 - (c1*a % f) - (c2*b % f)) / q
    
    # assert that the lattice vector v has the format
    # (u*q + c1*a % f + c2*b % f, a, b)
    # or, equivalently, that v = (v0, a, b) with
    # a*c1 + b*c2 = v0 in R/<f,q> (mod f and mod q)
    assert ((u*q + ((a*c1 + b*c2) % f)) == v0)

    return u, a, b


bound = 2
# polynomial with random degree and coefficients in [-bound, bound]
c1 = R.random_element(randint(0,N-1), -bound, bound) 
c2 = R.random_element(randint(0,N-1), -bound, bound) 

C1 = matrix_polynomial_product(c1)
C2 = matrix_polynomial_product(c2)

B = block_matrix(ZZ, [[q*Matrix.identity(N), C1, C2], [0, 1, 0], [0, 0, 1]])
print B
L = Lattice(B.transpose())

## Test if the lattice vectors have the expected format
for i in xrange(100):
    v = L.random_element() 
    u, a, b = lattice_point_to_polys(v)

print "Lattice points have the expected format."

print "c1 =", c1
print "c2 =", c2
v = L.shortest_vector()
print "Shortest (nonzero) vector of L:"
print v
u, a, b = lattice_point_to_polys(v)
print "a =", a
print "b =", b
print "u =", u

print "a*c1 + b*c2 % f =", ((a*c1 + b*c2) % f)
