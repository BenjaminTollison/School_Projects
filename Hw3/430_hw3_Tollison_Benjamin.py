### Benjamin Tollison Hw3 ###
import sympy as sp
from IPython.display import display, Math, Latex
from sympy.solvers.solveset import linsolve
from sympy import  lambdify, Matrix, sinh, cosh, exp, pi, symbols
import numpy as np
import matplotlib.pyplot as plt
def displayEquations(LHS,RHS):
    left = sp.latex(LHS)
    right = sp.latex(RHS)
    display(Math(left + '=' + right))
#### Analytical solution ####
# Governing equation
A,B,C,D,alpha_1,alpha_2,x,h,L = symbols('A B C D alpha_1 alpha_2 x h L')
u_1 = A*sinh(alpha_1*x) + B*cosh(alpha_1*x)
u_2 = C*sinh(alpha_2*x) + D*cosh(alpha_2*x)
P = 2*pi*x
bc1 = u_1.subs(x,0) - 100
bc2 = u_1.subs(x,1/2) - u_2.subs(x,0)
bc3 = u_2.diff(x).subs(x,1/2) - h*P*u_2.subs({x:1/2})
bc4 = u_1.diff(x).subs(x,1/2)  - h*P*(u_1.subs({x:1/2})-u_2.subs({x:0}))
dof = sp.nonlinsolve([bc1,bc2,bc3,bc4],[A,B,C,D])
(a_sol,b_sol,c_sol,d_sol) = next(iter(dof))
u_1 = u_1.subs({A:a_sol,B:b_sol})
u_2 = u_2.subs({C:c_sol,D:d_sol})
displayEquations('u_1',u_1)
displayEquations('u_2',u_2)
displayEquations('A',a_sol)
displayEquations('B',b_sol)
displayEquations('C',c_sol.subs({P:2*pi*x}))
displayEquations('D',d_sol.subs({P:2*pi*x}))
# Lets set alpha_1 = .2 and alpha_2 = .8, h = .5
from numpy import cosh,sinh,pi,exp
r = .1
u_1_graph = lambdify((x,alpha_1),u_1)
u_2_graph = lambdify((x, alpha_1, alpha_2, h), u_2)
x1_test_np_values = np.linspace(0,.5*r,10)
x2_test_np_values = np.linspace(.5*r,r,10)
plt.plot(x1_test_np_values,u_1_graph(x1_test_np_values,.2))
plt.plot(x2_test_np_values,u_2_graph(x1_test_np_values,.2,.8,.5))
plt.show()