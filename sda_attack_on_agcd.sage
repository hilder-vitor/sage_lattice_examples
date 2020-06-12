#	Trying to analyze the cost of the Simultaneous Diophantine approximation attack (SDA)
# against the AGCD problem.
#
#	For more details about this attack, we refer to the papers
# Algorithms for the Approximate Common Divisor Problem, by Galbraith, Gebregiyorgis, and Murphy;
# and  Fully homomorphic encryption over the integers, by Dijk, Gentry, Halevi, and Vaikuntanathan.
#



from sage.modules.free_module_integer import IntegerLattice as Lattice

def sample_r(rho):                                                                                                                                                                                           
	return ZZ.random_element(-2^rho, 2^rho+1)

def sample_q(gamma, eta):                                                                                                                                                                                    
	return ZZ.random_element(0, 2^(gamma - eta))

def sample_agcd(gamma, eta, rho, p):                                                                                                                                                                         
	return p * sample_q(gamma, eta) + sample_r(rho)

def create_basis_SDA(t, x0, x, rho):
	B = Matrix(ZZ, t+1, t+1); B[0, 0] = 2^(rho + 1);
	for i in range(1, t+1):
		B[0, i] = x[i-1]
		B[i, i] = -x0
	return B

def sample_x0(gamma, eta, rho, p):
	x0 = 2
	while(log(x0, 2) < gamma - 1):
		x0 = sample_agcd(gamma, eta, rho, p);

	return x0

def run(lam, k, t = 1):
	rho = lam
	eta = rho + k
	gamma = max(2*eta, ceil(lam * k^2 / log(lam, 2)))

	if 1 == t:
		t = ceil(2 * (gamma - rho) / (eta - rho))

	print("gamma = %d,  eta = %d,  rho = %d,  k = %d,  t = %d" % (gamma, eta, rho, k, t))

	p = random_prime(2^eta, lbound=2^(eta - 1))
	Zp = ZZ.quo(p)

	x0 = sample_x0(gamma, eta, rho, p)
	r0 = Zp(x0).lift_centered(); q0 = (x0 - r0) / p; assert(x0 == p*q0 + r0)

	## sample the AGCD instances used to construct the (t+1)-dimensional lattice
	x = vector([sample_agcd(gamma, eta, rho, p) for i in range(t)])
	r = vector([Zp(xi).lift_centered() for xi in x])
	q = (x - r) / p

	B = create_basis_SDA(t, x0, x, rho)
	L = Lattice(B)
	L.BKZ()
	det = L.volume()
	v = L.shortest_vector()
	lambda1 = v.norm()
	assert (lambda1 <= (sqrt(t+1) * (det ** (1/(t+1)))))

	target = vector([q0*2^(rho+1)] + [q0*r[i] - q[i]*r0 for i in range(t)])

	ratios = [(target.norm() / v.norm()).n()] + [(target.norm() / L.reduced_basis[i].norm()).n() for i in range(t+1)]
	print(ratios)

	mult = 1
	for i in range(t+1):
		vi = L.reduced_basis[i]
		if (target.norm() / vi.norm()) in ZZ:
			mult = vi[0] / target[0]
			print("target:        %s " % target)
			print("v%d:            %s " % (i, vi))
			print("mult:          %s " % mult)
			print("mult * target: %s " % vi) 

	print("")
	print("target == shortest vector? %s" % (target == v)) 
	print("target in BKZ base? %s" % (target in L.reduced_basis)) 
	print("mult * target in BKZ base? %s" % (mult * target in L.reduced_basis)) 


lam = 100 # security parameter
k = 3 # gap between p and noises
t = 10 # the dimention of the lattice is t+1
run(lam, k)
