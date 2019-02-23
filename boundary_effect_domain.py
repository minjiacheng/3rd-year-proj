import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
from functools import reduce
import random
from scipy import stats

N = 50 #define the dimension of the matrix
original_N = N
iteration = N #no. of iterations or columns for each simulation
#possible_resistors = [np.log(2), np.log(10), np.log(10)]
#below are param to generate random list of resistors
xmax = 5.0
xmin = -5.0
nstep = 1000
I = np.identity(N)
U_1 = 1000
trigger = 0

def f(x):
    if trigger == 0:
        return x**2
    elif trigger == 1:
        return 25-x**2
    
def prob_list_gen(xmin,xmax,nstep):
    x_val_list = []
    e_x_val_list = []
    y_val_list = []
    prob_list = []
    final_list = []
    
    for i in range(nstep):
        x = xmin + (xmax - xmin) * float(i) / nstep
        x_val_list = np.append(x_val_list,x)
        e_x_val_list = np.append(e_x_val_list,np.exp(x))
        y = f(x)
        y_val_list = np.append(y_val_list,y)
    for i in range(nstep):
        prob = y_val_list[i]/sum(y_val_list)
        #prob = 1/nstep
        #this turns it into a uniform distribution
        prob_list = np.append(prob_list,prob)
    #print(x_val_list)
    #print(prob_list)
    prob_list = np.round(prob_list,5)
    for i in range(len(prob_list)):
        prob_list[i]=prob_list[i]*100000
        temp_list = [x_val_list[i] for j in range(int(prob_list[i]))]
        final_list = np.append(final_list,temp_list)
    return final_list

possible_resistors = prob_list_gen(xmin,xmax,nstep)

def resistor_distribution():
    #define your random distribution here
    resistor = np.exp(random.choice(possible_resistors))
    return resistor

def random_generation():
    ##############################################################################
    #auto generate random V and H as initial configuration (v_i and h_i in range 0 to 10)
    v = np.zeros(N)
    i = 0
    while i < N:
        v[i] = resistor_distribution()
        i += 1
    V = np.zeros(shape=(N,N))
    V[0,0] = 1/v[0]
    i = 1
    while i < N:
        V[i,i] = 1/v[i] + 1/v[i-1]
        i += 1
    i = 0
    while i < N-1:
        V[i,i+1] = -1/v[i]
        V[i+1,i] = -1/v[i]
        i += 1
    H = np.zeros(shape=(N,N))
    i = 1
    while i < N:
        H[i,i] = resistor_distribution()
        i += 1
    h = np.zeros(N)
    j = 0
    while j < N:
        h[j] = H[j,j]
        j = j+1
    V = np.matrix(V)    
    H = np.matrix(H)
    ##############################################################################
    return v, h, V, H

def simulate():
    v, h, V, H = random_generation()
    V_list = [V]*N
    H_list = [H]*N
    V_resistors = [v]*N
    H_resistors = [h]*N
    A = V
    #print("A_1 is\n", A)
    i = 1
    invL = np.array([])
    sigma = np.array([])
    #all values of 1/L and sigma_N are stored in those arrays for plotting
    while i < iteration:
        v, h, V, H = random_generation()
        V_list[i] = V
        H_list[i] = H
        V_resistors[i] = v
        H_resistors[i] = h
        i += 1
        A = V + A * np.linalg.inv(I+ H*A)
        invA = np.linalg.inv(A)
        sigma_N = 1/invA[0,0]
        sigma_N = sigma_N / i
        invL = np.append(invL, [1/i])
        sigma = np.append(sigma, [sigma_N])
 
    #find resistance of network
    #print("sigma_N is \n", sigma_N)
    effective_R = sigma_N * N 
    effective_R = 1 / effective_R
    #for justification on this calc for effective resistance
    #see note book p14
    #print("overall_resistance is \n", overall_resistance)
    predicted_resistance = np.mean(possible_resistors)
    predicted_resistance = np.exp(predicted_resistance)
    #print(predicted_resistance,effective_R)
    percentage_deviation = ((effective_R-predicted_resistance)/effective_R)*100
    percentage_deviation = abs(percentage_deviation)
    #print("actual R:",effective_R)
    #print("predicted R:",predicted_resistance)
    return V_resistors, H_resistors, A, V_list, H_list, percentage_deviation

error_plot = []
xmins = []
while xmin < 0.0:
    list_of_percentage_error = []
    for i in range(1000):
        V_resistors, H_resistors, A, V_list, H_list, percentage_deviation = simulate()
        #find expected effective resistance
        #print("% deviation: ", percentage_deviation)
        list_of_percentage_error = np.append(list_of_percentage_error,percentage_deviation)
    avg_error = sum(list_of_percentage_error)/len(list_of_percentage_error)
    error_plot = np.append(error_plot,avg_error)
    xmins = np.append(xmins, xmin)
    xmin = xmin+0.5
    xmax = -xmin
    if xmin != 0.0:
        possible_resistors = prob_list_gen(xmin,xmax,nstep)
   
print("done half")
xmax = 5.0
xmin = -5.0
nstep = 1000
I = np.identity(N)
trigger = 1
possible_resistors = prob_list_gen(xmin,xmax,nstep)

error_plot2 = []
xmins2 = []
while xmin < 0.0:
    list_of_percentage_error = []
    for i in range(1000):
        V_resistors, H_resistors, A, V_list, H_list, percentage_deviation = simulate()
        #find expected effective resistance
        #print("% deviation: ", percentage_deviation)
        list_of_percentage_error = np.append(list_of_percentage_error,percentage_deviation)
    avg_error = sum(list_of_percentage_error)/len(list_of_percentage_error)
    error_plot2 = np.append(error_plot2,avg_error)
    xmins2 = np.append(xmins2, xmin)
    xmin = xmin+0.5
    xmax = -xmin
    if xmin != 0.0:
        possible_resistors = prob_list_gen(xmin,xmax,nstep)

plt.plot(abs(xmins), error_plot,label='f(x)=x^2')
plt.plot(abs(xmins2), error_plot2,label='f(x)=25-x^2')
plt.xlabel("bounds")
plt.ylabel("average % error over 1000 simulations")
plt.legend(loc='best')
plt.title("50 by 50 network. vary xmin and xmax")
