# -*- coding: utf-8 -*-
######################################################################################
#Replication of the G matrix computation in Section 5.2.1 using analytical derivatives
#in arbitrary precision using mpmath
#######################################################################################

import sympy
from sympy.matrices import zeros, eye

#Arbitrary precision using mpmath
import mpmath as mp
#from mpmath import mpf
#from mpmath import chop

prec=34

#prec=50  #uncomment to redo in 50-digit precision

mp.mp.dps = prec #set 34 or 50 digits of precision

#Set up parameters: Case (i): AMPF
symbol_str = 'alpha, phi_tau, phi_r, beta, gamma, sigma_r, sigma_tau'
alpha, phi_tau, phi_r, beta, gamma, sigma_r, sigma_tau = sympy.symbols(symbol_str, real=True)
ny, npara, ns = 2, 7, 4
Phi_1 = sympy.Matrix([[0, 0, -1/alpha*phi_r, 0],
[0, 1/beta - gamma*(1/beta - 1), 0, -phi_tau*(1/beta-1)],
[0, 0, 0, 0],
[0, 0, 0, 0]])
Phi_varepsilon = sympy.Matrix([[-1/alpha - 1/alpha**2*phi_r, 0],
[1/beta*(1/alpha+1/alpha**2*phi_r), -(1/beta-1)],
[1,0],
[0,1]])
Sigma = zeros(ny,ny)
Sigma[0,0] = sigma_r**2
Sigma[1,1] = sigma_tau**2

#Selection matrix A
A = zeros(ny,ns)
A[0,0] = 1
A[1,1] = 1


from IPython.display import Latex
sympy.init_printing(use_latex="mathjax")
ppi1, b1, e__r_t1, e__tau_t1 = sympy.symbols('\\hat{\pi}_{t-1}, \\hat{b}_{t-1}, e^r_{t-1}, e^\\tau_{t-1}')
ppi, b, e__r_t, e__tau_t = sympy.symbols('\\hat{\pi}_{t}, \\hat{b}_{t}, e^r_{t}, e^\\tau_{t}')
v1 = sympy.Matrix([[ppi1, b1, e__r_t1, e__tau_t1]]).transpose()
v0 = sympy.Matrix([[ppi, b, e__r_t, e__tau_t]]).transpose()
e0 = sympy.Matrix([[e__r_t, e__tau_t]]).transpose()
rhs = Phi_1 * v1 + Phi_varepsilon * e0
Latex(r'\begin{align*}\hat b_t &= ' + sympy.latex(rhs[1]) + r'\\ \hat\pi_t &= ' + sympy.latex(rhs[0])+'\end{align*}')

#compute symbolic spectral density
omega = sympy.symbols('omega', real=True)
x = sympy.exp(-sympy.I * omega)
H = A * (eye(ns) - Phi_1 * x).inv() * Phi_varepsilon
f = 1/(2*sympy.pi) * H * Sigma * H.H
f

para_list = [alpha, beta, gamma, phi_r, phi_tau, sigma_r, sigma_tau]
fp = sympy.Matrix(ny**2, npara, lambda i, j: f.vec()[i].diff(para_list[j]))
Z = (fp.H * fp)
Latex(r"\begin{dmath*} " + sympy.latex(Z[0,0]) + r"\end{dmath*}")



calibration = {alpha: mp.mpf(1.5),
               beta: mp.mpf(0.9804),
               gamma: mp.mpf(1.2),
               phi_r: mp.mpf(0.5),
               phi_tau: mp.mpf(0.5),
               sigma_tau: mp.mpf(1),
               sigma_r: mp.mpf(1)}

Z = Z.subs(calibration)

Gmp = mp.zeros(npara, npara)
for i in range(npara):
    for j in range(npara):
        fmpnumeric = sympy.lambdify(omega, Z[i,j], modules='mpmath')
        Gmp[i,j] = mp.quadgl(fmpnumeric, [-mp.pi, mp.pi],maxdegree=8) #NB: 768 quadrature nodes used here (m=8) due to rule of thumb of 3*(2^m) points

mp_eigs = mp.eig(Gmp)[0]
print('Mpmat quad  precision result')
print('-----------------------')
print('Eigenvalues of G: ', mp.chop(mp_eigs))