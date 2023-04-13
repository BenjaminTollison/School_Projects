### Benjamin Tollison ###

import numpy as np
import matplotlib.pyplot as plt

### problem 5 ###
E = 1.0 * 10**7
A = 0.1
EA = E * A
L = 10.0
g = 0.03
P = -10.0
# Creating the functions
def u_exact(x):
    u =  (1/EA) * (P * (x-L) + (1/20) * g * (L**5 - x**5)) + 2
    return u, 'u_exact'
def F_exact(x):
    return P - .25 * g * x**4, 'F_exact'
def u_approx(x):
    u = 2 + (P / EA) * (x-L) + (3 / (4 * EA)) * (((P * L + P) / L**2) + g * L**2 * ((L/6)-(1/4))) * (x**2 - L**2)
    return u, 'u_approx'
def F_approx(x):
    F = P + (3/2) * x * (((P*L + P) / L**2) + g * L**2 *((L/6) - (1/4)))
    return F, 'F_approx'
def u_variationally_consistant(x):
    u =  (1 / EA) * (P * x - (1/20) * g * x**5) + 2
    return u, 'u_vc'
def F_variationally_consistant(x):
    return P - (1/4) * g * x**4, 'F_vc'
#creating a function to plot all the graphs to compare them
def plot(function1, function2):
    x_values = np.linspace(0, L, 100)
    plt.title('{} vs. {}'.format(function1(1)[1], function2(1)[1]))
    plt.plot(x_values, function1(x_values)[0],label = '{}'.format(function1(1)[1]))
    plt.plot(x_values, function2(x_values)[0],label = '{}'.format(function2(1)[1]))
    plt.xlabel('x')
    plt.legend()
    plt.ylabel('{}(y) and {}(y)'.format(function1(1)[1], function2(1)[1]))
    plt.show()
#comparing the exact to fem
plot(u_approx,u_exact)
plot(F_approx,F_exact)
#comparing Exact to Variationally_Consistant
plot(u_variationally_consistant,u_exact)
plot(F_variationally_consistant,F_exact)
### Problem 6 ###
L = 10
E = 10**7
A = 0.1
EA = E * A
def F_exact(x):
    return (200/L**3) * (.25 * (x**4 +L**4) - x -L) , 'F_exact'
def u_exact(x):
    return (10/(EA * L**3)) * (x**5 - 10 * x**2), 'u_exact'
def u_approx(x):
    return (200/(EA*L**3)) * (-2990*x + (199/2)* x**2), 'u_approx'
def F_approx(x):
    return (200/(L**3)) * (-2990 + (199)* x), 'F_approx'
def F_variationally_consistant(x):
    F = (200 / L**3) * (.25*x**4 - x) + 970/3
    return F, 'F_variationally_consistant'
def u_variationally_consistant(x):
    u = (200 / (EA * L**3)) * ((1/20) * (x**5) - .5 * x**2) + (970 / (3 * EA)) * x
    return u, 'u_variationally_consistant'
plot(u_approx,u_exact)
plot(F_approx,F_exact)
plot(u_variationally_consistant,u_exact)
plot(F_variationally_consistant, F_exact)
### Problem 7
L = 10
EA = 1*10**6
def u_exact(x):
    return (100 / (EA * L**3)) * ((1/10)*x**5 - x**2 + (L - (1/10)*L**4)*x), 'u_exact'
def F_exact(x):
    return (100 / L**3) * (((1/2)*x**4 - 2*x) + (L-(1/10)*L**4)), 'F_exact'
def u_approx(x):
    return (1/EA) * (284.4642 * x - 882.839*x**2 + 65.489*x**3), 'u_approx'
def F_approx(x):
    return 284.464 - 2*882.839*x + 65.489*3*x**2, 'F_approx'
def u_variationally_consistant(x):
    return (200/(EA*L**3)) * ((1/20)*((x**5)) - .5*(x**2)) + (99/EA)*x , 'u_vc'
def F_variationally_consistant(x):
    return (200/L**3)*(.25*x**4 - x ) +99, 'F_vc'
plot(u_approx,u_exact)
plot(F_approx,F_exact)
plot(u_variationally_consistant,u_exact)
plot(F_variationally_consistant, F_exact)