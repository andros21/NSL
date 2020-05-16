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
print('$ Psi2 (cpp) $')
print_ccode(GS_norm2, standard='C99')

print()
print('$ Psi2 (python) $')
print(str(GS_norm2).replace('exp', 'np.exp').replace(
    'pi', 'np.pi').replace('sqrt', 'np.sqrt'))

print()
print('$ Integ (cpp) $')
print_ccode(simplify(H_GS/GS), standard='C99')
