# -*- coding: utf-8 -*-
######################################################################################
#Replication of the G matrix computation in Section 5.2.1 using analytical derivatives
#in double precision
######################################################################################
#

import numpy as np
import sympy
from sympy.matrices import zeros, eye

#Precision digits
prec1=34
prec2 = 50

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

calibration = {alpha: sympy.Float('1.5', precision=prec1),
               beta: sympy.Float('0.9804', precision=prec1),
               gamma: sympy.Float('1.2', precision=prec1),
               phi_r: sympy.Float('0.5', precision=prec1),
               phi_tau: sympy.Float('0.5',precision=prec1),
               sigma_tau: sympy.Float('1', precision=prec1),
               sigma_r: sympy.Float('1', precision=prec1)}

#Substitute parameter and use Gaussian quadrature with 500 nodes
Z = Z.subs(calibration)
x, w = np.polynomial.legendre.leggauss(500)
t = 0.5*(x + 1)*(np.pi*2) - np.pi # integral on [-pi,pi]
G = np.zeros((npara, npara))
#Converting the weights to [-pi,pi]
w=w*np.pi
for i in range(npara):
    for j in range(npara):
        fnumeric = sympy.lambdify(omega, Z[i,j])
        G[i,j] = (w * np.real(fnumeric(t))).sum() 
np.set_printoptions(precision=3)

#Eigenvalues of G
np_eigs = np.linalg.eigvals(G)
print('Numpy double precision result')
print('-----------------------')
print('Eigenvalues of G: ', np_eigs)


