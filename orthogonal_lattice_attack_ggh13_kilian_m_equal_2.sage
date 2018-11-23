from sage.modules.free_module_integer import IntegerLattice as Lattice

def q_ary_orthogonal_lattice(A, q):
    '''
    Return the lattice given by the integer vectors that are orthogonal to A
    over Zq, that is, 
    L(A)_q_ort := {u in Z^n : A*u = 0 mod q}
    '''
    m = A.nrows()
    n = A.ncols()
    if m >= n:
        return False

    Zq = IntegerModRing(q)
    A1 = Matrix(Zq, A[:, 0:m])
    invA1 = ((A1)^-1).lift() # assuming that the m fisrt collums of A are LI over Zq
    A2 = A[:, m:n]

    B = block_matrix([[q, -invA1 * A2], [0, 1]])

    return Lattice(B.transpose(), lll_reduce = True) # B_ort is already an LLL-reduced basis


# return a polynomial h in Zx such that h*g = 1 over Rq,
# or, equivalently,  ((h*g mod f) mod q) = 1 over ZZ['x']
def inverse_polynomial(g):
    return Zx((1 / Rq(g)).lift().coefficients(sparse=False))


def phi(h): # isomorphism from Zx/f to Z^N (return a vector representing h)
    v = h.coefficients(sparse=False)
    return v + [0]*(N - len(v))

def inv_phi(vec_h): # isomorphism from Z^N to ZZ['x']/f (return the polynomial represented by vec_h)
    return Zx([ai for ai in vec_h])

# Return a matrix H such that for any p in Zx/<f>, H*phi(p) equals h*p % f,
# that is, the matrix that represents the product by the polynomial h.
def matrix_polynomial_product(h):
    H = Matrix(ZZ, [phi(h * x^i % f) for i in xrange(N)])
    return H.transpose()

def centered_mod(h):
    vec_h = h.coefficients(sparse=False)
    return Zx([Zq(ai).lift_centered() for ai in vec_h])

def size(poly):
    return vector(ZZ, phi(poly)).norm(p=Infinity)

# All the vectors in L are supposed to have the format 
#     phi(a1) || phi(a2)|| ... || phi(a_M+1)
# for M+1 polynomials a1, ..., a_M+1 of Zx/<f>.
# This function receives a v from L and return such polynomials.
def lattice_point_to_polys(v):
    return [inv_phi(v[i*N:(i+1)*N]) for i in xrange(M+1)]
    
def encode(m):
    r = Zx.random_element(randint(N/2,N-1), -B_E, B_E)
    return centered_mod(((r*g + m) * inv_z) % f)

###################################################################################
############################    MAIN PART       ###################################

# Use M < N < q
M = 2
N = 2**3  # Degree of polynomial used as modulus and dimension of lattice
q = random_prime(N^2 * M, proof=None, lbound=N) # arbitrary upper bound
Zx.<x> = ZZ['x']
f = x^N + 1     # modulus defining Zx / <f>
R = Zx.quotient(f)

B_g = N
B_M = N
B_E = N

g = Zx.random_element(randint(N/2,N-1), -B_g, B_g)
Zq = GF(q)
Rq = Zq['x'].quotient(f)

print "Ring: ZZ[x]/<", f, ">"
print "Generator of the principal ideal: g =", g
print "q =", q


z = Zx.random_element(randint(0,N-1), int(-q/2), int(q/2))
inv_z = inverse_polynomial(z)
assert 1 == ((z * inv_z % f) %q)


K = Matrix.random(Zq, M, M).lift()
assert(0 != (K.det() % q)) # K is invertible over Zq

# vector of polynomials with random degree and coefficients in [-B_M, B_M]
m1 = [Zx.random_element(randint(0,g.degree()-1), -B_M, B_M) for i in xrange(M)]
c1 = [encode(m1_i) for m1_i in m1]
c1K = (vector(c1) * K % q) # randomizing the encodings
m2 = [Zx.random_element(randint(0,g.degree()-1), -B_M, B_M) for i in xrange(M)]
c2 = [encode(m2_i) for m2_i in m2]
c2K = (vector(c2) * K % q)
m3 = [Zx.random_element(randint(0,g.degree()-1), -B_M, B_M) for i in xrange(M)]
c3 = [encode(m3_i) for m3_i in m3]
c3K = (vector(c3) * K % q)

print ""
print "Vectors of encodings:"
print "c1*K =", c1K 
print "c2*K =", c2K 
print "c3*K =", c3K 

C11 = matrix_polynomial_product(c1K[0])
C12 = matrix_polynomial_product(c1K[1])
C21 = matrix_polynomial_product(c2K[0])
C22 = matrix_polynomial_product(c2K[1])
C31 = matrix_polynomial_product(c3K[0])
C32 = matrix_polynomial_product(c3K[1])

A = block_matrix(ZZ, [[C11, C21, C31], [C12, C22, C32]])

L = q_ary_orthogonal_lattice(A, q)

a,b,c = lattice_point_to_polys(L.reduced_basis[0])
print ""
print "Orthogonal lattice L created."
print "(aprox) shortest: (a, b, c) =", (a, b, c)
print "size(a, b, c) =", max(size(a), size(b), size(c))

for i in xrange(100):
    v = L.random_element()
    a, b, c = lattice_point_to_polys(v)

    # assert that the polynomials "in L" are orthogonal to the original encodings (without Kilian's randomization)
    assert(0 == Rq(a*c1[0] + b*c2[0] + c*c3[0]))
    assert(0 == Rq(a*c1[1] + b*c2[1] + c*c3[1]))

print "Vectors in L are orthogonal to the vectors of encodings [c11, c21, c31] and [c12, c22, c32]"

