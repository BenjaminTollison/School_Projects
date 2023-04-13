import sympy as sp
from IPython.display import display, Math, Latex
from sympy.solvers.solveset import linsolve
from sympy import  lambdify, Matrix 
import numpy as np
import matplotlib.pyplot as plt
def displayEquations(LHS,RHS):
    left = sp.latex(LHS)
    right = sp.latex(RHS)
    display(Math(left + '=' + right))
#### Problem 9 ####
    # Exact solution
x,L,EA,Q,a,b,alpha,Delta_T = sp.symbols('x L EA Q a b alpha Delta_T')
D2u = -1/EA * 200*(1-x**3/L**3) + alpha * Delta_T
Du = sp.integrate(D2u,x) + a
u = sp.integrate(Du,x) + b
fbc = EA * Du.subs(x,L) - Q
kbc = u.subs(x,0)
dof = sp.linsolve([fbc,kbc],[a,b])
(a_sol,b_sol) = next(iter(dof))
u = u.subs({a:a_sol,b:b_sol})
displayEquations('u(x)',u)
strain = u.diff(x)
displayEquations('strain(x)',strain)
#Graphing the exact solution
u_np = u.subs({Delta_T:400,alpha:10**-6,EA:10**6,L:10,Q:-500})
u_np = lambdify(x,u_np)
x_linspace = np.linspace(0,10,1000)
plt.plot(x_linspace,u_np(x_linspace))
plt.xlabel('x')
plt.ylabel('u(x)')
plt.title('Displacement Exact solution')
plt.show()
# strain
strain_np = strain.subs({Delta_T:400,alpha:10**-6,EA:10**6,L:10,Q:-500})
strain_np = lambdify(x,strain_np)
plt.plot(x_linspace,strain_np(x_linspace))
plt.xlabel('x')
plt.ylabel('epison(x)')
plt.title('Strain Exact solution')
plt.show()
# FEM soltution
    # displaying equations to make handwork easier
def N_i(h):
    return sp.Matrix([[1-x/h],[x/h]])
def F_T(h):
    return sp.integrate(EA * alpha * Delta_T * N_i(h).diff(x), (x, 0, h))
def F_dist(h):
    f = 200*(1-x**3/L**3)
    return sp.integrate(f * N_i(h),(x,0,h))
def K(h):
    return EA/h * sp.Matrix([[1,-1],[-1,1]])
    # 0-1 sympy
print('Displaying matrices for nodes 0 to 1')
displayEquations('K_0', K(5))
displayEquations('N_0', N_i(5))
displayEquations('F_T_0', F_T(5))
displayEquations('F_dist_0', F_dist(5))
LHS_0 = F_T(5) + F_dist(5)
print('printing LHS_0')
displayEquations('LHS_0', LHS_0)
LHS_1 = F_T(3) + F_dist(3)
displayEquations('LHS_1', LHS_1)
LHS_2 = F_T(2) + F_dist(2) + sp.Matrix([[0],[-500]])
displayEquations('LHS_2', LHS_2)
# Assemblying the full LHS
LHS = sp.Matrix([LHS_0,LHS_1, LHS_2])
displayEquations('LHS', LHS)
# plugging in the values
LHS = LHS.subs({EA:10**6,Delta_T:500,alpha:10**-6,L:10})
displayEquations('LHS', LHS)
# Direct Assembly 
K_0 = K(5).subs({EA:10**6})
K_1 = K(3).subs({EA:10**6})
K_2 = K(2).subs({EA:10**6})
K_global = np.zeros((4,4))
K_global[0][0],K_global[0][1] = K_global[1][0] + K_0[0], K_global[0][1] + K_0[1]
K_global[1][0],K_global[1][1] = K_global[1][0]+K_0[2], K_global[1][1]+K_0[3]
K_global[1][1],K_global[1][2] = K_global[1][1] + K_1[0], K_global[1][2] + K_1[1]
K_global[2][1],K_global[2][2] = K_global[2][1]+K_1[2], K_global[2][2]+K_1[3]
K_global[2][2],K_global[2][3] = K_global[2][2] + K_2[0], K_global[2][3] + K_2[1]
K_global[3][2],K_global[3][3] = K_global[3][2]+K_2[2], K_global[3][3]+K_2[3]
print(K_global)
LHS_direct = np.array([[LHS[0]], [LHS[1]+LHS[2]], [LHS[3]+LHS[4]],[LHS[5]]])
print(LHS_direct)
# applying kinematic bc of u(0) = 0
K_global[0][0:] = 0
K_global[1][0] = 0
for i in range(0,3):
    K_global[i][i] = 1
print(K_global)
LHS_direct[0] = 0
displacement_array = np.linalg.inv(K_global) @ LHS_direct
print(displacement_array)
# comparing fem u(x) to analytical solution
x_postions = np.array([0,5,8,10])
plt.plot(x_postions,displacement_array,label='FEM solution')
plt.plot(x_linspace,u_np(x_linspace),label='analytical solution')
plt.xlabel('x')
plt.ylabel('u(x)')
plt.title('Comparing solutions')
plt.legend()
plt.show()
y_prime_list = []
for i in range(1, len(displacement_array)):
    y_prime_list.append((displacement_array[i] - displacement_array[i-1])/(x_postions[i] - x_postions[i-1]))
print(y_prime_list)
y0 = [y_prime_list[0]] * 10
y1 = [y_prime_list[1]] * 10
y2 = [y_prime_list[2]] * 10
x0 = np.linspace(0, 5, 10)
x1 = np.linspace(5, 8, 10)
x2 = np.linspace(8, 10, 10)
plt.plot(x0, y0)
plt.plot(x1, y1)
plt.plot(x2, y2)
plt.plot(x_linspace,strain_np(x_linspace))
plt.xlabel('x')
plt.ylabel('epison(x)')
plt.title('Strain solutions')
plt.show()