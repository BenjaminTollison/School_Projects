import numpy as np
import matplotlib.pyplot as plt


def f(var):
    x = var
    A = 17/12
    f = x**4 / 12 - x**2 + A * x
    return f

i = np.arange(0,3,0.01)
di = 0.01

plt.plot(i, f(i),'r')
plt.hlines(0,0,3,'k')
plt.xlabel('x')
plt.ylabel('f(x)')
plt.title('Analytical Solution')
plt.show()

h = 0.2

xp = []
j = np.arange(0,3+h,h)
for z in j:
    k = z
    k = round(k,2)
    xp.append(k)

f_prime = []
f_double_prime = []
for i in range(len(xp)):
    f_prime.append((2/3) * xp[i])
    f_double_prime.append(((2/3) + xp[i]))

plt.plot()

### FDM ###
    # case a
# building the coefficient matrix for different mesh sizes
n = np.arange(2,10,1)
delta_x = lambda nodes: 1 / (nodes - 1)
alpha = 6
kap = 2 + alpha**2*delta_x**2
for i in n:
    Coefficient_Matrix = np.zeros((i,i))
    Coefficient_Matrix[i-1][i-2:i] = [-1, kap]
    Coefficient_Matrix[0][0:2] = [kap, -1]
    for i in range(1,i-1):
        Coefficient_Matrix[i][i-1:i+2] = [-1, kap, -1]
def u_fdm(number_of_nodes):
    i = number_of_nodes
    u_0 = 5
    u_1 = 100
    Coefficient_Matrix = np.zeros((i,i))
    Coefficient_Matrix[i-1][i-2:i] = [-1, kap]
    Coefficient_Matrix[0][0:2] = [kap, -1]
    for j in range(1,i-1):
        Coefficient_Matrix[j][j-1:j+2] = [-1, kap, -1]
    # building the Boundary conditions array
    boundary_conditions = np.zeros(i)
    boundary_conditions[0] = u_0
    boundary_conditions[i-1] = u_1
    # finding the nodes inbetween the boundaries
    u_array = np.array([0,4,14,38,100])
    return np.matmul(Coefficient_Matrix,u_array), Coefficient_Matrix
print(u_fdm(5))