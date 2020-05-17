from sympy import Symbol, integrate, oo, Derivative, simplify, exp
from sympy.printing import print_ccode

x = Symbol('x')
mu = Symbol('mu', real=True, positive=True)
sigma = Symbol('sigma', real=True, positive=True)

Gp = exp(-(x-mu)**2/(2*sigma)**2)
Gm = exp(-(x+mu)**2/(2*sigma)**2)
V = x**4 - (5/2) * x**2

GS = Gp + Gm

GS_norm2 = integrate(GS**2, (x, -oo, oo))

H_GS = -(1/2) * Derivative(GS, x, 2).doit() + V * GS

GS_norm2 = simplify(GS**2/GS_norm2)

print()
print('##### For main.cpp code and 08.ipynb #####', end='\n\n')

print('$ Psi2 (cpp) $')
print_ccode(GS_norm2, standard='C99')

print()
print('$ Psi2 (python) $')
print(str(GS_norm2).replace('exp', 'np.exp').replace(
    'pi', 'np.pi').replace('sqrt', 'np.sqrt'))

print()
print('$ Integ (cpp) $')
print_ccode(simplify(H_GS/GS), standard='C99')

print()
print('##### For QMC_1D Code #####', end='\n\n')

val = Symbol('val')
v = Symbol('v')
Gp = exp(-(v-mu)**2/(2*sigma)**2)
Gm = exp(-(v+mu)**2/(2*sigma)**2)
GS = Gp + Gm
V = val**4 - (5/2) * val**2

print('$ Potential 0st der (cpp) $')
print_ccode(simplify(V), standard='C99')

print()
print('$ Potential 1st der (cpp) $')
print_ccode(simplify(Derivative(V, val, 1).doit()), standard='C99')

print()
print('$ Potential 2st der (cpp) $')
print_ccode(simplify(Derivative(V, val, 2).doit()), standard='C99')

print()
print('$ psi 0st der (cpp) $')
print_ccode(simplify(GS), standard='C99')

print()
print('$ psi 2st der (cpp) $')
print_ccode(simplify(Derivative(GS, v, 2).doit()), standard='C99')
