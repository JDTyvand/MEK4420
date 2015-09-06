from numpy import *
from matplotlib.pyplot import *

a = 1
b = 1
num_segments = 100
num_points = num_segments+1
rad_pos = linspace(0, 2*pi, num_points)
x = a*cos(rad_pos)
y = b*sin(rad_pos)
A_array = []
B_array = []

for i in range(num_segments):
	A_array.append(array([x[i], y[i]]))
	B_array.append(array([x[i+1], y[i+1]]))

class Segment:

	def __init__(self, A, B):
		self.A = A
		self.B = B
		self.midpoint = (self.A+self.B)/2
		self.dx = self.B[0] - self.A[0]
		self.dy = self.B[1] - self.A[1]
		self.dl = linalg.norm(self.B - self.A)
		self.nx = -self.dy/self.dl
		self.ny = self.dx/self.dl
		self.n = [self.nx, self.ny]
		self.dpdn = self.nx + self.ny
		self.lhs = zeros(num_segments)
		self.rhs = 0
		self.r1 = []
		self.r1_norm = zeros(num_segments)
		self.r2 = []
		self.r2_norm = zeros(num_segments)

	def calculate(self, segments):
		for i in range(len(segments)):
			if segments[i] is self:
				self.lhs[i] = -pi
				self.r1.append(0)
				self.r2.append(0)
			else:
				self.r1.append(segments[i].A - self.midpoint)
				self.r2.append(segments[i].B - self.midpoint)
				self.r1_norm[i] = linalg.norm(self.r1[i])
				self.r2_norm[i] = linalg.norm(self.r2[i])
				self.lhs[i] = -arccos(dot(self.r1[i],self.r2[i])/(self.r1_norm[i]*self.r2_norm[i]))
				self.rhs = self.rhs + segments[i].dpdn*((log(self.r1_norm[i]) +log(self.r2_norm[i]))*segments[i].dl)/2


segments = [Segment(A_array[i],B_array[i]) for i in range(num_segments)]

for i in range(num_segments):
	segments[i].calculate(segments)

A = zeros((num_segments,num_segments))
for i in range(num_segments):
	A[i] = segments[i].lhs

B = [segments[i].rhs for i in range(num_segments)]
phi = linalg.solve(A,B)

print -(a*a*x[10])/(x[10]**2 + y[10]**2)
print phi[10]


