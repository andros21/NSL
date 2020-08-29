#!/usr/bin/env python
# -*- coding: utf-8 -*-

# THEORETICAL VALUES

import numpy as np


def Vpot(x):
    return (x**2 - 2.5)*x**2
    # return 0.5*x**2


hbar = 1
m = 1
a = 6
N = 1000  # number of iterations
# Step sizes
x = np.linspace(-a/2, a/2, N)
dx = x[1] - x[0]  # the step size
V = Vpot(x)
# The central differences method: f" = (f_1 - 2*f_0 + f_-1)/dx^2
CDiff = np.diag(np.ones(N-1), -1)-2*np.diag(np.ones(N), 0) + \
                np.diag(np.ones(N-1), 1)
# np.diag(np.array,k) construct a "diagonal" matrix using the np.array
# The default is k=0. Use k>0 for diagonals above the main diagonal,
# and k<0 for diagonals below the main diagonal
# Hamiltonian matrix
H = (-(hbar**2)*CDiff)/(2*m*dx**2) + np.diag(V)
# Compute eigenvectors and their eigenvalues
E, psi = np.linalg.eigh(H)
# Take the transpose & normalize
psi = np.transpose(psi)
psi = psi/np.sqrt(dx)
