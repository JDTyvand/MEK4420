from numpy import *
from matplotlib.pyplot import *
import time
import warnings
warnings.filterwarnings('ignore')
set_printoptions(precision=4, suppress=True)



def ellipse(a,b):
	"""
	Function calculating the added mass 
	coefficients for an ellipse,
	or a circle if a=b.
	"""
	print
	if a == b:
		print('Case: circle with radius %d' % a)
	else:
		print('Case: ellipse with a = %d, b = %d' % (a,b))
	def calculate(i):
		#Distance from midpoint of segment i to 
		#start point of each of the other segments
		r1 = linalg.norm(array([x[:-1],y[:-1]]).T - 
				array([midpointx[i],midpointy[i]]),axis=1)   
		#Distance from midpoint of segment i to 
		#end point of each of the other segments
		r2 = linalg.norm(array([x[1:],y[1:]]).T - 
				array([midpointx[i],midpointy[i]]),axis=1)
		#Array of angles between the distance vectors 
		#to the other segments  
		theta = -arccos((dl**2 - r2**2 - r1**2)/(-2*r2*r1))	     
		#Calculates the right-hand side integral 
		#for x, y, and rotation
		rhs11 = sum(n1*(log(r1)+log(r2))*0.5*dl)	             
		rhs22 = sum(n2*(log(r1)+log(r2))*0.5*dl)		     
		rhs66 = sum(n6*(log(r1)+log(r2))*0.5*dl)	
		#Adds the angles to the matrix A	    
		A[i] = theta	
		#Adds rhs to the B-arrays					    
		B11[i] = rhs11										 
		B22[i] = rhs22										 
		B66[i] = rhs66										 

	start = time.time()	
	#Number of segments									
	N = 4000			
	#Array of N+1 points alomng the geometry									 
	rad_pos = linspace(0, 2*pi, N+1)	
	#X and Y-values of the points					
	x = a*cos(rad_pos)										
	y = b*sin(rad_pos)			
	# X and Y-values of the midpoints							
	midpointx = (x[1:] + x[:-1])/2							 
	midpointy = (y[1:] + y[:-1])/2		
	#Length of each segment					 
	dl = linalg.norm(array([x[1:],y[1:]]) - 				 
			array([x[:-1],y[:-1]]), axis=0)
	#Components of the normal vectors of segments
	n1 = -(y[1:] - y[:-1])/dl 						       	 		 
	n2 = (x[1:] - x[:-1])/dl 								               
	n6 = (midpointx*n2 - midpointy*n1)                      

	#Matrix and arrays for storing lhs and rhs of equations
	A = zeros((N,N))                                         
	B11 = zeros(N)                                          
	B22 = zeros(N)											
	B66 = zeros(N)											 

	for i in range(N):
		calculate(i)

	#Fills the diagonal of A to remove potential NaN values
	fill_diagonal(A,-pi)									

	#Calculates phi for the three directions
	phi11 = linalg.solve(A,B11)								 
	phi22 = linalg.solve(A,B22)								 
	phi66 = linalg.solve(A,B66)                              

	#Calculates the added mass coefficients
	m11 = sum(phi11*n1*dl)                                   
	m22 = sum(phi22*n2*dl)									 
	m66 = sum(phi66*n6*dl)									 

	#Assembles a matrix of added mass coefficients for printing
	M = zeros((3,3))										 
	M[0][0] = m11
	M[1][1] = m22
	M[2][2] = m66

	print M

	end = time.time()
	print ('Calculated in %4f seconds' % (end-start))

	if a == b:
		"""
		This loop calculates the exact 
		solution at the midpoints of 
		segments, and compares with 
		the numerical solution. 
		This is also plotted.
		"""
		exact = -(a**2*midpointx)/(midpointx**2 
					+ midpointy**2)
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
	This function is very similar 
	to the one for the ellipse, 
	and comments are only
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
		#Adds -pi to the index for the current segment 
		theta[i] = -pi 	
		#Changes rest of NaN values to 0									 
		theta[isnan(theta)] = 0								 
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
	#Ensures segments is divisable by 4
	N = N/4 * 4 											 
	Np = N/4 												 
	x = zeros(N+1)
	y = zeros(N+1)
	for i in range(Np+1):
		"""
		This loops adds points to all 
		four sides at the same time,
		using cosine spacing to get 
		denser spacing in the corners.
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

	dl = linalg.norm(array([x[1:],y[1:]]) - 
		  array([x[:-1],y[:-1]]), axis=0)
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