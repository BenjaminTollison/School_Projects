### Intialize Main ###
def FEA_prototype(Text_file):
    import numpy as np
    import sympy as sp
    import csv
    f = open(Text_file)
    reader = csv.reader(f, delimiter=' ')
    data_string = []
    for row in reader:
        data_string.append(row) # in a string form
    # read Mesh
    node_num = []
    x_distance = []
    for i in range(1, int(data_string[0][0]) + 1):
        node_num.append(int(data_string[i][0])-1)
        x_distance.append(float(data_string[i][1]))
    # read properties
    element_num = []
    connectivity_of_element = []
    for i in range(int(data_string[0][0]) + 1 , int(data_string[0][0]) + int(data_string[0][1]) + 1):
        element_num.append(float(data_string[i][0]))
        connectivity_of_element.append((int(data_string[i][1])-1,int(data_string[i][2])-1))
    for i in range(int(data_string[0][0]) + int(data_string[0][1]) + 1, int(data_string[0][0]) + int(data_string[0][1]) + 2):
        E = float(data_string[i][0])
        height = float(data_string[i][1])
        width = float(data_string[i][2])
    # read constraints
    num_of_constraints = int(data_string[int(data_string[0][0]) + int(data_string[0][1]) + 2][0])
    dof_to_be_zeros = []
    for i in range(num_of_constraints):
        dof_to_be_zeros.append(int(data_string[int(data_string[0][0]) + int(data_string[0][1]) + 3][i])-1)
    # read loads
    number_of_point_loads = int(data_string[len(data_string) - 2][0])
    dof = int(data_string[len(data_string) - 1][0])
    load = float(data_string[len(data_string) - 1][1])
    # get Element DOF
    def get_element_dof(based_on_node_number_list):
        q = map(chr, range(0,(based_on_node_number_list)))
        a = []
        for i in bytearray(range(97,97 + (based_on_node_number_list))).decode("utf-8"):
            a += i
        a = sp.symbols(a)
        a = sp.Matrix(a)
        for i in dof_to_be_zeros:
            a[i-1] = 0
        return a
    # get Element K
    def get_element_k(nodes_that_are_connnected):
        EA = E * height * width
        x,y = nodes_that_are_connnected
        h = abs(x_distance[(x)-1]-x_distance[(y)-1]) # difference in x
        K = EA / h * np.array([[1, -1], [-1, 1]])
        return K
    # assemble Global Stiffness Matrix with direct assembly
        # "Busted like yo gf, scream something back" -GySgt Maley
    global_stiffness_matrix = np.zeros((len(node_num), len(node_num)))
    for i in range(0, len(connectivity_of_element)):
        x,y = connectivity_of_element[i]
        global_stiffness_matrix[x-1][x-2], global_stiffness_matrix[x-1][x-1] =  get_element_k(connectivity_of_element[i])[0]
        global_stiffness_matrix[y-1][y-2], global_stiffness_matrix[y-1][y-1] =  get_element_k(connectivity_of_element[i])[1]
    output_information = [node_num,x_distance,element_num,connectivity_of_element,E,height,width,num_of_constraints,dof_to_be_zeros,number_of_point_loads,dof,load]
        # only returning Part 1 in python numbering system
    output_text_file = open("Output_{}.txt".format(Text_file), "w")
    output_text_file.write("#################################\n"+
                           "Benjamin Tollison FEA output file\n"+
                           "#################################\n"+
                           "Number of nodes used: {}\n".format(len(node_num))+
                           "Node list: {}\n".format(node_num)+
                           "x position of nodes: {}\n".format(x_distance)+
                           "Number of Elements: {}\n".format(len(element_num))+
                           "Connection of Elements: {}\n".format(connectivity_of_element)+
                           "Young's Modulus, height, and width of rod: {}, {}, {}\n".format(E,height,width)+
                           "Number of constraints: {}\n".format(num_of_constraints)+
                           "Degrees of freedom to be zero: {}\n".format(dof_to_be_zeros)+
                           "Number of point loads: {}\n".format(number_of_point_loads)+
                           "DOF: {}\n".format(dof)+
                           "Load: {}\n".format(load))
    output_text_file.close()
    f.close()
    return None
# Executing the code from the selected file 
FEA_prototype('data_1.txt')