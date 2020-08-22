#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Theoretical Curves

import numpy as np

T = np.linspace(0.5, 2.0, num=100)
beta = 1/T
h = 0.02
b = beta
J = 1.0
Ns = 50
th = np.tanh(J/T)
thN = th**Ns
ch = 1/th
U = -J*(th + ch*thN)/(1 + thN)
X = beta*np.exp(2*beta*J)*(1-thN)/(1+thN)
C = ((beta*J)**2)*(((1+thN+(Ns-1)*(th**2)+(Ns-1)*(ch**2)*thN)/(1+thN)) -
                   Ns*((th+ch*thN)/(1+thN))**2)
l1 = np.exp(b*J)*np.cosh(b*h)+np.sqrt(np.exp(2*b*J) *
                                      np.cosh(b*h)*np.cosh(b*h)-2*np.sinh(2*b*J))
l2 = np.exp(b*J)*np.cosh(b*h)-np.sqrt(np.exp(2*b*J) *
                                      np.cosh(b*h)*np.cosh(b*h)-2*np.sinh(2*b*J))
Z = l1**Ns + l2**Ns
M = (np.exp(b*J)*np.sinh(b*h)*((l1**(Ns-1))*(1+np.exp(b*J)*np.cosh(b*h)/np.sqrt(np.exp(2*b*J)*np.cosh(b*h)*np.cosh(b*h)-2 *
                                                                                np.sinh(2*b*J)))+(l2**(Ns-1))*(1-np.exp(b*J)*np.cosh(b*h)/np.sqrt(np.exp(2*b*J)*np.cosh(b*h)*np.cosh(b*h)-2*np.sinh(2*b*J)))))/(Z)
