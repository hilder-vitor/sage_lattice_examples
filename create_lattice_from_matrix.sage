from sage.modules.free_module_integer import IntegerLattice as Lattice

n = 4

B = Matrix.random(ZZ, n, n)

print "Basis of L:"
print B

# create a lattice with basis given by the rows of B
L = Lattice(B, lll_reduce=False) # default: lll_reduce = True

for i in xrange(B.nrows()):
    row_i = vector(B[i, :])
    assert (row_i in L) == True

print "All rows of B belong to L"
print ""

for j in xrange(B.ncols()):
    col_j = vector(B[:, j])
    if col_j not in L:
        print "Column %d of B doesn't belong to L:" % (j+1)
        print col_j
        break

