from numpy import *
from math import ceil
from matplotlib.pyplot import *

N = 100

d1 = -1
d2 = 1
x = zeros(N)
t = linspace(d1, d2, 100)
for i in range(N):
	x[i] = d1 + (d2-d1)/2*(1-cos(pi/(N-1)*i))
plot(t, x, 'ro')
show()