"""
Compute the number of moulins our domain should have using
observation-constrained density
"""

import numpy as np
import scipy.stats

# Parameters for surface elevation equation
a = 100e3
b = 5e3
k1 = 6
z0 = 390
surf = lambda x: k1*(np.sqrt(x+b) - np.sqrt(b)) + z0

# Empirically derived moulin density
pdf = lambda x: 42.12111317*scipy.stats.norm.pdf(x, loc=1138.25114462, scale=280.12484981)

# parameters for integration
dy = 25e3
dx = 1e3
x = np.arange(0, 100e3, dx)
z = surf(x)

n_moulins = 0
for (i, xi) in enumerate(x):
    zmean = z[i]
    # Number moulins = density (#/km2) * area (km2)
    n_moulins += pdf(zmean)*dx*dy/1e6

print(n_moulins)
