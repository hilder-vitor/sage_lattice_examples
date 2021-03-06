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
M = 5
N = 2**3  # Degree of polynomial used as modulus and dimension of lattice
q = random_prime(N^2 * M, proof=None, lbound=N) # arbitrary upper bound

assert(M < N < q)

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

print "M =", M

z = Zx.random_element(randint(0,N-1), int(-q/2), int(q/2))
inv_z = inverse_polynomial(z)
assert 1 == ((z * inv_z % f) %q)


K = Matrix.random(Zq, M, M).lift()
assert(0 != (K.det() % q)) # K is invertible over Zq

# vector of polynomials with random degree and coefficients in [-B_M, B_M]
m = [[Zx.random_element(randint(0,g.degree()-1), -B_M, B_M) for j in xrange(M)] for i in xrange(M+1)]
c = [[encode(m[i][j]) for j in xrange(M)] for i in xrange(M+1)]
cK = [(vector(c[i]) * K % q) for i in xrange(M+1)]  # randomizing the encodings

print ""
print "Vectors of encodings:"
for i in xrange(M+1):
    print "c%i*K = %s" %(i, cK[i])

C = [[matrix_polynomial_product(cK[i][j]) for i in xrange(M+1)] for j in xrange(M)]

A = block_matrix(ZZ, C)#[[C11, C21, C31], [C12, C22, C32]])
L = q_ary_orthogonal_lattice(A, q)
assert(L.rank() == (M+1) * N)     # is the rank equal to N*(M+1)?
assert(L.rank() == L.dimension()) # is the lattice full-rank?
assert(L.volume() == q^(M*N))     # is the determinant equal to q^(M*N)?

print ""
print "Orthogonal lattice L created."
print "Rank = dimension = N*(M+1) = ", L.rank()

for i in xrange(150):
    v = L.random_element()
    polys = lattice_point_to_polys(v)

    # assert that the "polynomials in L" are orthogonal to the original encodings (without Kilian's randomization)
    for j in xrange(M):
        inner_prod = 0
        for i in xrange(M+1):
            inner_prod += polys[i] * c[i][j]
        assert(0 == Rq(inner_prod))

print "Vectors in L are orthogonal to the vectors of encodings."

