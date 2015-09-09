from numpy import *
from matplotlib.pyplot import *
import time
import warnings
warnings.filterwarnings('ignore')
set_printoptions(precision=4, suppress=True)



def ellipse(a,b):

	def calculate(i):
		r1 = linalg.norm(array([x[:-1],y[:-1]]).T - array([midpointx[i],midpointy[i]]),axis=1)
		r2 = linalg.norm(array([x[1:],y[1:]]).T - array([midpointx[i],midpointy[i]]),axis=1)
		theta = -arccos((dl**2 - r2**2 - r1**2)/(-2*r2*r1))
		rhs11 = sum(nx*(log(r1)+log(r2))*0.5*dl)
		rhs22 = sum(ny*(log(r1)+log(r2))*0.5*dl)
		rhs66 = sum(n66*(log(r1)+log(r2))*0.5*dl)
		A[i] = theta
		B11[i] = rhs11
		B22[i] = rhs22
		B66[i] = rhs66

	start = time.time()
	a = a
	b = b
	num_segments = 4000
	num_points = num_segments+1
	rad_pos = linspace(0, 2*pi, num_points)
	x = a*cos(rad_pos)
	y = b*sin(rad_pos)
	midpointx = (x[1:] + x[:-1])/2
	midpointy = (y[1:] + y[:-1])/2
	dl = linalg.norm(array([x[1:],y[1:]]) - array([x[:-1],y[:-1]]), axis=0)
	nx = -(y[1:] - y[:-1])/dl
	ny = (x[1:] - x[:-1])/dl
	r = array([midpointx,midpointy]).T
	n = array([nx, ny]).T
	n66 = cross(r,n)

	A = zeros((num_segments,num_segments))
	B11 = zeros(num_segments)
	B22 = zeros(num_segments)
	B66 = zeros(num_segments)

	for i in range(num_segments):
		calculate(i)

	fill_diagonal(A,-pi)

	phi11 = linalg.solve(A,B11)
	phi22 = linalg.solve(A,B22)
	phi66 = linalg.solve(A,B66)

	m11 = sum(phi11*nx*dl)
	m22 = sum(phi22*ny*dl)
	m66 = sum(phi66*n66*dl)

	M = zeros((3,3))
	M[0][0] = m11
	M[1][1] = m22
	M[2][2] = m66

	print M

	end = time.time()
	print ('Calculated in %4f seconds' % (end-start))

def square(a):

	def calculate(i):
		r1 = linalg.norm(array([x[:-1],y[:-1]]).T - array([midpointx[i],midpointy[i]]),axis=1)
		r2 = linalg.norm(array([x[1:],y[1:]]).T - array([midpointx[i],midpointy[i]]),axis=1)
		theta = -arccos((dl**2 - r2**2 - r1**2)/(-2*r2*r1))
		theta[i] = -pi
		theta[isnan(theta)] = 0
 		rhs11 = sum(nx*(log(r1)+log(r2))*0.5*dl)
		rhs22 = sum(ny*(log(r1)+log(r2))*0.5*dl)
		rhs66 = sum(n66*(log(r1)+log(r2))*0.5*dl)
		A[i] = theta
		B11[i] = rhs11
		B22[i] = rhs22
		B66[i] = rhs66

	start = time.time()
	d1 = -a
	d2 = a
	N = 4000
	N = N/4 * 4
	Np = N/4
	x = zeros(N+1)
	y = zeros(N+1)
	for i in range(Np+1):
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
	nx = -(y[1:] - y[:-1])/dl
	ny = (x[1:] - x[:-1])/dl
	r = array([midpointx,midpointy]).T
	n = array([nx, ny]).T
	n66 = cross(r,n)

	A = zeros((N,N))
	B11 = zeros(N)
	B22 = zeros(N)
	B66 = zeros(N)

	for i in range(N):
		calculate(i)

	phi11 = linalg.solve(A,B11)
	phi22 = linalg.solve(A,B22)
	phi66 = linalg.solve(A,B66)

	m11 = sum(phi11*nx*dl)
	m22 = sum(phi22*ny*dl)
	m66 = sum(phi66*n66*dl)

	M = zeros((3,3))
	M[0][0] = m11
	M[1][1] = m22
	M[2][2] = m66

	print M

	end = time.time()
	print ('Calculated in %4f seconds' % (end-start))

if __name__ == '__main__':
    ellipse(2,2)
    ellipse(2,1)
    square(1)