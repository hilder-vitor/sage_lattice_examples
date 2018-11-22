from sage.modules.free_module_integer import IntegerLattice as Lattice

######################################################################
#   This script includes methods to generate the q-ary lattices
#       L(A)_q = { u in Z^n : u = A^T v mod q for some v in Z^m}
#       L(A)_q_ort = { u in Z^n : 0 = Au mod q }
# for any given matrix A in Z^mxn, with m < n (under the assumption that
# the first m collums of A are linearly independent).   
#######################################################################


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


def q_ary_lattice(A, q):
    '''
    Return the lattice given by the mod q integer linear combinations of the rows
    of A, that is,
    L(A)_q := {u in Z^n : there exists v in Z^m : u = v*A mod q}
    '''
    m = A.nrows()
    n = A.ncols()
    if m >= n:
        return False

    Zq = IntegerModRing(q)
    A1 = Matrix(Zq, A[:, 0:m])
    invA1 = ((A1)^-1).lift() # assuming that the m fisrt collums of A are LI over Zq
    A2 = A[:, m:n]

    B = block_matrix([[1, 0], [(invA1 * A2).transpose(), q]])

    return Lattice(B.transpose(), lll_reduce = True) # B_ort is already an LLL-reduced b

# Use m < n < q
m = 4
n = 6
q = random_prime(10*n*m, proof=None, lbound=n) # arbitrary upper bound 10*n*m
Zq = IntegerModRing(q)

A = Matrix([[ZZ.random_element(-q//2, q//2) for i in xrange(n)] for i in xrange(m)])
print "A:"
print A
print "q =", q

Lq_ort = q_ary_orthogonal_lattice(A, q)

### Asserting that vectors in Lq_orthogoal are really orthogonal to the collums of A over Zq
for i in xrange(n):
    bi = Lq_ort.reduced_basis[i] # we just need to test for the basis vectors
    assert(0 == A*bi % q)

assert(Lq_ort.volume() == q**m)
print "det(Lq_ort) =", Lq_ort.volume(), "is q^m, as expected."


Lq = q_ary_lattice(A, q)
### Asserting that for every vector u of Lq there exist an integer vector v such that u = A^T * v mod q
for i in xrange(n):
    bi = Lq.reduced_basis[i] # we just need to test for the basis vectors
    v = Matrix(Zq, A.transpose()).solve_right(bi)
    assert(bi % q == A.transpose()*v % q)

assert(Lq.volume() == q**(n-m))
print "det(Lq) =", Lq.volume(), "is q^(n - m), as expected."
