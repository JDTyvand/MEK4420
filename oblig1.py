from numpy import *
from matplotlib.pyplot import *
import time
import warnings
warnings.filterwarnings('ignore')
set_printoptions(precision=4, suppress=True)



def ellipse(a,b):
	"""
	Function calculating the added mass coefficients for an ellipse,
	or a circle if a=b.
	"""
	print
	if a == b:
		print('Case: circle with radius %d' % a)
	else:
		print('Case: ellipse with a = %d, b = %d' % (a,b))
	def calculate(i):
		r1 = linalg.norm(array([x[:-1],y[:-1]]).T - 
				array([midpointx[i],midpointy[i]]),axis=1)   #Distance from midpoint of segment i to start point of each of the other segments
		r2 = linalg.norm(array([x[1:],y[1:]]).T - 
				array([midpointx[i],midpointy[i]]),axis=1)   #Distance from midpoint of segment i to end point of each of the other segments
		theta = -arccos((dl**2 - r2**2 - r1**2)/(-2*r2*r1))	 #Array of angles between the distance vectors to the other segments
		rhs11 = sum(n1*(log(r1)+log(r2))*0.5*dl)			 #Calculates the right-hand side integral for x-direction
		rhs22 = sum(n2*(log(r1)+log(r2))*0.5*dl)			 #Calculates the right-hand side integral for y-direction
		rhs66 = sum(n6*(log(r1)+log(r2))*0.5*dl)			 #Calculates the right-hand side integral for rotation
		A[i] = theta										 #Adds the angles to the matrix A
		B11[i] = rhs11										 #Adds rhs for x-direction to array B11
		B22[i] = rhs22										 #Adds rhs for y-direction to array B22
		B66[i] = rhs66										 #Adds rhs for rotation to array B66

	start = time.time()										 #Starts counter for timing of calculations
	N = 4000												 #Number of segments
	rad_pos = linspace(0, 2*pi, N+1)						 #Array of N+1 points alomng the geometry
	x = a*cos(rad_pos)										 #X-values of the points
	y = b*sin(rad_pos)										 #Y-values of the points
	midpointx = (x[1:] + x[:-1])/2							 #X-values of the midpoints of the segments
	midpointy = (y[1:] + y[:-1])/2							 #Y-values of the midpoints of the segments
	dl = linalg.norm(array([x[1:],y[1:]]) - 				 #Length of each segment
			array([x[:-1],y[:-1]]), axis=0)
	n1 = -(y[1:] - y[:-1])/dl 						       	 #X-component of normal vector on segment		 
	n2 = (x[1:] - x[:-1])/dl 								 #Y-component of normal vector on segment              
	n6 = (midpointx*n2 - midpointy*n1)                       #Rotation component of normal vector on segment

	A = zeros((N,N))                                         #N-by-N matrix for storing of angle values
	B11 = zeros(N)                                           #Array for storing rhs-values in x-direction
	B22 = zeros(N)											 #Array for storing rhs-values in y-direction
	B66 = zeros(N)											 #Array for storing rhs-values in rotation direction

	for i in range(N):
		calculate(i)

	fill_diagonal(A,-pi)									 #Fills the diagonal of A to remove potential NaN values

	phi11 = linalg.solve(A,B11)								 #Calculates values of phi in the x-direction
	phi22 = linalg.solve(A,B22)								 #Calculates the values of phi in the y-direction
	phi66 = linalg.solve(A,B66)                              #Calculates the values of phi in rotation direction

	m11 = sum(phi11*n1*dl)                                   #Calculates the added mass coefficient m11
	m22 = sum(phi22*n2*dl)									 #Calculates the added mass coefficient m22
	m66 = sum(phi66*n6*dl)									 #Calculates the added mass coefficient m66

	M = zeros((3,3))										 #Assembles a matrix of added mass coefficients for printing
	M[0][0] = m11
	M[1][1] = m22
	M[2][2] = m66

	print M

	end = time.time()
	print ('Calculated in %4f seconds' % (end-start))

	if a == b:
		"""
		This loop calculates the exact solution at the midpoints of segments, and compares with 
		the numerical solution. This is also plotted.
		"""
		exact = -(a**2*midpointx)/(midpointx**2 + midpointy**2)
		err = abs(exact-phi11).max()
		print('The maximum error between exact and numerical solution is %.5f' % err)
		plot(exact, 'b-', label='exact solution')
		plot(phi11, 'r-', label='numerical solution')
		xlabel('segment')
		ylabel('phi')
		legend(loc='upper right', numpoints = 1)
		savefig('exact_vs_numerical.png')

def square(a):
	"""
	This function is very similar to the one for the ellipse, and comments are only
	given for the differing functionality
	"""
	print
	print('Case: square with sides %d' %(2*a))
	def calculate(i):
		r1 = linalg.norm(array([x[:-1],y[:-1]]).T - 
				array([midpointx[i],midpointy[i]]),axis=1)
		r2 = linalg.norm(array([x[1:],y[1:]]).T - 
				array([midpointx[i],midpointy[i]]),axis=1)
		theta = -arccos((dl**2 - r2**2 - r1**2)/(-2*r2*r1))
		theta[i] = -pi 										 #Adds -pi to the index for the current segment to compensate for NaN-values
		theta[isnan(theta)] = 0								 #Changes NaN values to 0
 		rhs11 = sum(n1*(log(r1)+log(r2))*0.5*dl)
		rhs22 = sum(n2*(log(r1)+log(r2))*0.5*dl)
		rhs66 = sum(n6*(log(r1)+log(r2))*0.5*dl)
		A[i] = theta
		B11[i] = rhs11
		B22[i] = rhs22
		B66[i] = rhs66

	start = time.time()
	d1 = -a
	d2 = a
	N = 4000
	N = N/4 * 4 											 #Integer division an multiplication to ensure segments is divisable by 4
	Np = N/4 												 
	x = zeros(N+1)
	y = zeros(N+1)
	for i in range(Np+1):
		"""
		This loops adds points to all four sides at the same time,
		using cosine spacing to get denser spacing in the corners.
		"""
		x[i] = d1 + (d2-d1)/2*(1-cos(pi/(Np)*i))
		y[i] = -a
		x[i+Np] = a
		y[i+Np] = d1 + (d2-d1)/2*(1-cos(pi/(Np)*i))
		x[i+2*Np] = -(d1 + (d2-d1)/2*(1-cos(pi/(Np)*i)))
		y[i+2*Np] = a
		x[i+3*Np] = -a
		y[i+3*Np] = -(d1 + (d2-d1)/2*(1-cos(pi/(Np)*i)))

	midpointx = (x[1:] + x[:-1])/2
	midpointy = (y[1:] + y[:-1])/2

	dl = linalg.norm(array([x[1:],y[1:]]) - array([x[:-1],y[:-1]]), axis=0)
	n1 = -(y[1:] - y[:-1])/dl
	n2 = (x[1:] - x[:-1])/dl
	n3 = zeros_like(n1)
	n6 = (midpointx*n2 - midpointy*n1)   

	A = zeros((N,N))
	B11 = zeros(N)
	B22 = zeros(N)
	B66 = zeros(N)

	for i in range(N):
		calculate(i)

	phi11 = linalg.solve(A,B11)
	phi22 = linalg.solve(A,B22)
	phi66 = linalg.solve(A,B66)

	m11 = sum(phi11*n1*dl)
	m22 = sum(phi22*n2*dl)
	m66 = sum(phi66*n6*dl)

	M = zeros((3,3))
	M[0][0] = m11
	M[1][1] = m22
	M[2][2] = m66

	print M

	end = time.time()
	print ('Calculated in %4f seconds' % (end-start))

if __name__ == '__main__':
    """
    Main function running three different cases
    """
    ellipse(1,1)
    ellipse(2,1)
    square(1)