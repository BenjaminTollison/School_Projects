#### Benjamin Tollison ####

# I created my own code because I don't feel that if I just
# changed the alpha values in the other codes that I would 
# actually learn the material. 

import numpy as np
import matplotlib.pyplot as plt

### Analytical Solutions ###
    # case a
x = np.linspace(0,1,1000)
alpha = np.array([0.25, 0.50,0.73, 1., 2., 4.,6.33, 8., 10.73])
u_a = lambda alpha, x: (100 / np.sinh(alpha)) * np.sinh(alpha * x)
    #plotting u(x)
for i in alpha:
    plt.plot(x, u_a(i,x), label='' + str(i) + ' \u03B1')
plt.title('Case A exact Temperature')
plt.xlabel('x')
plt.ylabel('u(x) in Celcius')
plt.legend()
plt.show()
    # Q(a)
kap = .02
r = 0.1
A = np.pi * r**2
Q_a = lambda alpha, x: -kap * A * alpha * (100 / np.sinh(alpha)) * np.cosh(alpha * x)
for i in alpha:
    plt.plot(x, Q_a(i,x), label='' + str(i) + ' \u03B1')
plt.title('Case A exact Bulk heat loss')
plt.xlabel('x')
plt.ylabel('Q(x)')
plt.legend()
plt.show()
    # case b
    #u(b)
u_b = lambda alpha, x: (100 / np.cosh(alpha)) * np.cosh(alpha * x)
for i in alpha:
    plt.plot(x, u_b(i,x), label='' + str(i) + ' \u03B1')
plt.title('Case B exact Temperature')
plt.xlabel('x')
plt.ylabel('u(x) in Celcius')
plt.legend()
plt.show()
    #Q(b)
Q_b = lambda alpha, x: -kap * A * alpha * (100 / np.cosh(alpha)) * np.sinh(alpha * x)
for i in alpha:
    plt.plot(x, Q_a(i,x), label='' + str(i) + ' \u03B1')
plt.title('Case B exact Bulk heat loss')
plt.xlabel('x')
plt.ylabel('Q(x)')
plt.legend()
plt.show()

    # Case c
def u_c(alpha,x):
    r = 0.1
    h = 1
    P = np.pi * r**2
    u_x = (50 * h * P * np.exp(-alpha)) / (alpha * np.cosh(alpha)**2) * np.sinh(alpha * x) + 100 / np.cosh(alpha) * np.cosh(alpha * x)
    return u_x

def Q_c(alpha, x):
    r = 0.1
    h = 1
    P = np.pi * r**2
    k = 1
    A = P
    u_prime = (50 * h * P * np.exp(-alpha)) / (np.cosh(alpha)**2) * np.cosh(alpha * x) - (100 * alpha) / np.cosh(alpha) * np.sinh(alpha * x)
    Q = k*A * u_prime
    return Q


#### graphing
for i in alpha:
    plt.plot(x, u_c(i,x), label='' + str(i) + ' \u03B1')
plt.title('Case C exact Temperature')
plt.xlabel('x')
plt.ylabel('u(x) in Celcius')
plt.legend()
plt.show()

for i in alpha:
    plt.plot(x, Q_c(i,x), label='' + str(i) + ' \u03B1')
plt.title('Case C Bulk heat loss')
plt.xlabel('x')
plt.ylabel('Q(x)')
plt.legend()
plt.show()

    # creating table values for exact value graphs
u_a_table = np.zeros([len(alpha), 11])
o=-1
for i in alpha:
     o += 1
     p = 0
     for j in np.arange(0,1.1,.1):
        u_a_table[o,p] = '%.1f' % u_a(i,j)
        p += 1
print(u_a_table)

u_a_table = np.zeros([len(alpha), 11])
o=-1
for i in alpha:
     o += 1
     p = 0
     for j in np.arange(0,1.1,.1):
        u_a_table[o,p] = '%.1f' % u_a(i,j)
        p += 1
print(u_a_table)

u_c_table = np.zeros([len(alpha), 11])
o=-1
for i in alpha:
     o += 1
     p = 0
     for j in np.arange(0,1.1,.1):
        u_c_table[o,p] = '%.1f' % u_c(i,j)
        p += 1
print(u_c_table)

# Making Coefficent matrix for case a
def Case_a_Coeff_matrix(number_of_nodes,alpha):
    i = number_of_nodes
    delta_x = 1 / (i+1)
    kap = 2 + alpha**2 * delta_x**2
    Coefficient_Matrix = np.zeros((i,i))
    Coefficient_Matrix[i-1][i-2:i] = [-1, kap]
    Coefficient_Matrix[0][0:2] = [kap, -1]
    for j in range(1,i-1):
        Coefficient_Matrix[j][j-1:j+2] = [-1, kap, -1]
    return Coefficient_Matrix
    
# building the Boundary conditions array 
def Boundary_Conditions_array(number_of_nodes):
    i = number_of_nodes
    u_0 = 0
    u_1 = 100
    boundary_conditions = np.zeros(i)
    boundary_conditions[0] = u_0
    boundary_conditions[i-1] = u_1
    return (boundary_conditions)

def u_array(nodes,alpha):
    n = nodes
    u = np.dot(np.linalg.inv(Case_a_Coeff_matrix(n,alpha)), Boundary_Conditions_array(n))
    return u
# in order to get the full line you must account for the boundary conditions that are not 
# accounted for in the the u_array. The u_array only finds the values inbetween

def full_list_of_temperature_values(total_nodes,alpha):
    n = total_nodes - 2
    a = alpha
    temp_list_y = []
    temp_list_y.append(0)
    for i in u_array(n,a):
        temp_list_y.append(i)
    temp_list_y.append(100)
    return temp_list_y

alpha_list = [0.25, 0.50,0.73, 1., 2., 4.,6.33, 8., 10.73]
for a in alpha_list:
    alpha = a
    total_nodes = 23
    dx = 1 / (total_nodes-1)
    x_values_list = []
    for i in np.arange(0,1+dx,dx):
        x_values_list.append(i)
    plt.plot(x_values_list,full_list_of_temperature_values(total_nodes,a),label='' + str(a) + ' \u03B1')
    plt.title('Case A FDM solution')
    plt.ylabel('Celcius')
    plt.xlabel('x')
    plt.legend()
# Extracting the table values for case a fdm
# reminder the vaules outputed are the in-between values and you have to add 0/100 to the ends of the final table values
u_a_table = np.zeros([len(alpha_list), 11])
o=-1
for i in alpha_list:
     o += 1
     p = 0
     for j in np.arange(0,9,1):
        u_a_table[o,p] = '%.1f' % u_array(9,i)[j]
        p += 1
print('Case A fdm')
print(u_a_table)

### Making a general code for the FDM solution

def Coeff_b_matrix(number_of_nodes,alpha):
    i = number_of_nodes
    delta_x = 1 / (i+1)
    h = 1
    r = 0.1
    P = np.pi * r**2
    kap = 2 + alpha**2 * delta_x**2
    k_prime = kap / 2
    Coefficient_Matrix = np.zeros((i,i))
    Coefficient_Matrix[i-1][i-2:i] = [-1, kap]
    Coefficient_Matrix[0][0:2] = [k_prime, -1]
    for j in range(1,i-1):
        Coefficient_Matrix[j][j-1:j+2] = [-1, kap, -1]
    return Coefficient_Matrix

# building the Boundary conditions array 
def Boundary_Conditions_array(number_of_nodes):
    i = number_of_nodes
    u_prime_0 = 0
    u_1 = 100
    boundary_conditions = np.zeros(i)
    boundary_conditions[0] = u_prime_0
    boundary_conditions[i-1] = u_1
    return (boundary_conditions)
# re-defined u_array
def u_b_array(nodes,alpha):
    n = nodes
    u = np.dot(np.linalg.inv(Coeff_b_matrix(n,alpha)), Boundary_Conditions_array(n))
    return u
# changed the value of u_array for u_b_array
def full_list_of_temperature_values(total_nodes,alpha):
    n = total_nodes - 2
    a = alpha
    temp_list_y = []
    temp_list_y.append(0)
    for i in u_b_array(n,a):
        temp_list_y.append(i)
    temp_list_y.append(100)
    return temp_list_y
# same plotting method but with replaced labeling
# omitting the first (x,y) to get rid of edge effects caused by the prime function
for a in alpha_list:
    alpha = a
    total_nodes = 75
    dx = 1 / (total_nodes-1)
    x_values_list = []
    for i in np.arange(0,1+dx,dx):
        x_values_list.append(i)
    plt.plot(x_values_list[1:],full_list_of_temperature_values(total_nodes,a)[1:],label='' + str(a) + ' \u03B1')
    plt.title('Case B FDM solution')
    plt.ylabel('Celcius')
    plt.xlabel('x')
    plt.legend()
# Generating the table
# only add 100 to the end of table values
u_b_table = np.zeros([len(alpha_list), 11])
o=-1
for i in alpha_list:
     o += 1
     p = 0
     for j in np.arange(0,10,1):
        u_b_table[o,p] = '%.1f' % u_b_array(10,i)[j]
        p += 1
print('Case B fdm')
print(u_b_table)


### Case 3 FDM

def Coeff_c_matrix(number_of_nodes,alpha):
    i = number_of_nodes
    delta_x = 1 / (i+1)
    h = 1
    r = 0.1
    P = np.pi * r**2
    kap = 2 + alpha**2 * delta_x**2
    k_prime = ((h) * delta_x * P) + kap / 2
    Coefficient_Matrix = np.zeros((i,i))
    Coefficient_Matrix[i-1][i-2:i] = [-1, kap]
    Coefficient_Matrix[0][0:2] = [k_prime, -1]
    for j in range(1,i-1):
        Coefficient_Matrix[j][j-1:j+2] = [-1, kap, -1]
    return Coefficient_Matrix
# building the Boundary conditions array 
def Boundary_Conditions_array(number_of_nodes):
    i = number_of_nodes
    u_prime_0 = 0
    u_1 = 100
    boundary_conditions = np.zeros(i)
    boundary_conditions[0] = u_prime_0
    boundary_conditions[i-1] = u_1
    return (boundary_conditions)
# re-defined u_array
def u_c_array(nodes,alpha):
    n = nodes
    u = np.dot(np.linalg.inv(Coeff_c_matrix(n,alpha)), Boundary_Conditions_array(n))
    return u
# changed the value of u_array for u_c_array
def full_list_of_temperature_values(total_nodes,alpha):
    n = total_nodes - 2
    a = alpha
    temp_list_y = []
    temp_list_y.append(0)
    for i in u_c_array(n,a):
        temp_list_y.append(i)
    temp_list_y.append(100)
    return temp_list_y
# same plotting method but with replaced labeling
for a in alpha_list:
    alpha = a
    total_nodes = 500
    dx = 1 / (total_nodes-1)
    x_values_list = []
    for i in np.arange(0,1+dx,dx):
        x_values_list.append(i)
    plt.plot(x_values_list[1:],full_list_of_temperature_values(total_nodes,a)[1:],label='' + str(a) + ' \u03B1')
    plt.title('Case C FDM solution')
    plt.ylabel('Celcius')
    plt.xlabel('x')
    plt.legend()
# printing out the table
# only add 100 to end of values
u_c_table = np.zeros([len(alpha_list), 11])
o=-1
for i in alpha_list:
     o += 1
     p = 0
     for j in np.arange(0,10,1):
        u_c_table[o,p] = '%.1f' % u_c_array(10,i)[j]
        p += 1
print('Case C fdm')
print(u_c_table)

### calculating the bulk heat flux at Q_fdm(1) with taylor series expansion

# adjusting temperature function to take in account each case 
def full_list_of_temperature_values(total_nodes,alpha,case_array):
    n = total_nodes
    a = alpha
    temp_list_y = []
    temp_list_y.append(0)
    for i in case_array(n,a):
        temp_list_y.append(i)
    temp_list_y.append(100)
    return temp_list_y

    # finding bulk heat loss at the begining of the heat fin
def Q_fdm(case_array,nodes,alpha):
    u_n = full_list_of_temperature_values(nodes,alpha,case_array)[nodes-1]
    u_n_1 = full_list_of_temperature_values(nodes,alpha,case_array)[nodes-2]
    k = .02
    r = 0.1
    A = np.pi * r**2
    delta_x = 1 / (nodes)
    return - k * A * ((u_n - u_n_1) / delta_x) 

# finding the total heat disapating
nodes_num = [4,5,8,15,25,35,50,75,85,161]
for j in [u_array,u_b_array,u_c_array]:
    print(j)
    for i in nodes_num:
        print( round( 1 / (i),4), "%.5f" % (Q_fdm(j,i,alpha_list[4])))
# finding total heat from exact equations
def Q_exact(case_array,nodes,alpha):
    k = 1
    r = 0.1
    A = np.pi * r**2
    delta_x = 1 / (nodes)
    if case_array == u_array:
        return Q_a(alpha, 1)
    elif case_array == u_b_array:
        return Q_b(alpha, 1)
    else:
        return Q_c(alpha, 1)
    
for j in [u_array,u_b_array,u_c_array]:
    print(j)
    print( "%.5f" % (Q_exact(j,i,alpha_list[4])))

### Finding the rate of convergence, starting with the case 1 fdm and alpha = 2 ###
def beta(case_array):
    nodes_list = []
    b = []
    percent_error = []
    delta_x = []
    for i in range(1,13):
        nodes_list.append(2**i)
    for j in range(0,len(nodes_list)-1):
        b.append(round(np.log2((Q_exact(case_array,nodes_list[j],alpha_list[4]) - Q_fdm(case_array,nodes_list[j],alpha_list[4]))/(Q_exact(case_array,nodes_list[j+1],alpha_list[4]) - Q_fdm(case_array,nodes_list[j+1],alpha_list[4]))),5))
        percent_error.append(round(abs((Q_exact(case_array,nodes_list[j],alpha_list[4]) - Q_fdm(case_array,nodes_list[j],alpha_list[4])) / (Q_exact(case_array,nodes_list[j],alpha_list[4]))) * 100, 5))
        delta_x.append(1 / nodes_list[j])
    return b,percent_error, delta_x
print(beta(u_array))
print(beta(u_b_array))
print(beta(u_c_array))

#### finding percent error ###
def percent_error(case_array):
    nodes_list = []
    p_err = []
    for i in range(1,13): # reduce this 13 to something else if you need faster load times
        nodes_list.append(2**i)
    for j in range(0,len(nodes_list)-1):
        p_err.append(abs((Q_exact(case_array,nodes_list[j],alpha_list[4]) - Q_fdm(case_array,nodes_list[j],alpha_list[4])) / (Q_exact(case_array,nodes_list[j],alpha_list[4]))) * 100)
    return p_err
# These lists take awhile to calculate
print(percent_error(u_array))
print(percent_error(u_b_array))
print(percent_error(u_c_array))

### double checking analytical solutions with a symbolic solver
import sympy as sp
alpha, x, gamma = sp.symbols('alpha, x, gamma')
u = sp.symbols('u',cls=sp.Function)
u_0 = (100 - (100 / sp.cosh(alpha)) * sp.sinh(alpha)) / ((2 / gamma) * sp.cosh(alpha))
u = (100 / sp.cosh(alpha)) * sp.sinh(alpha * x) + ((2 * u_0) / gamma) * sp.cosh(alpha * x)
print(u.diff(x,1))

# finding u_c_exact via cramer's rule
J = sp.Matrix([[sp.exp(alpha), sp.exp(alpha)],[gamma, - gamma]])
C = sp.det(J)
C_1 = sp.det(sp.Matrix([[100, u_0],[sp.exp(- alpha), - gamma]]))
C_2 = sp.det(sp.Matrix([[sp.exp(alpha), gamma],[100, u_0]]))
A = C_1 / C
B = C_2 / C
u = A * sp.exp(alpha* x) + B * sp.exp(- alpha * x)
print(sp.simplify(u))
print(u.diff(x,1))