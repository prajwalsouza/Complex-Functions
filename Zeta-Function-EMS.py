from mpmath import mp
from scipy import special

# The following formula is based on Euler-Maclaurin Summation as mentioned in 
# https://math.dartmouth.edu/archive/m56s13/public_html/Nguyen_proj.pdf 
# which is the work of Hanh Nguyen. 

def zetaEMS(s, N, v):
	sum1 = mp.mpc(0)
	s = mp.mpc(s)
	for n in range(1, N):
		sum1 = sum1 + (n**(-s))
	
	sum1 = sum1 + ((N**(1 - s))/mp.mpc(s - 1))
	sum1 = sum1 + ((N**(-s))/mp.mpc(2))

	sum2 = 0
	for k1 in range(1, v + 1):
		t1 = (bernoulli(2*k1)/mp.mpc(fact(2*k1)))
		prd = mp.mpc(1)
		for h in range(0, (2*k1) - 1):
			prd = prd*(s + mp.mpc(h))
		t2 = prd
		t3 = N**(1 - s - (2*k1))
		sum2 = sum2 + (t1*t2*t3)

	return sum1 + sum2

# Bernoulli numbers calculated using Scipy's special function module

def bernoulli(val):
	return special.bernoulli(val)[val]

# Computing the factorial of a number

def fact(d):
	pd = 1
	for a in range(1, d + 1):
		pd = pd*a

	return pd

# Computing Combination nCr (Binomial Coefficients)

def nCr(n, r):
	return mp.mpf(fact(n))/mp.mpf(fact(r)*fact(n - r))


complexinput = 0.5+14.13j
output = zetaEMS(complexinput, N=100, v=100)
# As mentioned in the work, N has to chosen to be of the order of magnitude of the complex input. But, here, I've chosen a large value.
# value of v has been chosen depending on accuracy. (More than v = 100, bernoulli numbers are extremely high and can lead to errors)

print("Input : ")
print(complexinput)
print("Zeta Function : ")
print(output)
