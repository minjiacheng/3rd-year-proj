import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import skewnorm
from mpl_toolkits import mplot3d
from scipy.interpolate import interp1d
import random

def find_k(skewness):
    k = 4.0 / (skewness**2)
    return k

N = 5 #define the dimension of the matrix
iteration = N #no. of iterations or columns for each simulation
I = np.identity(N)
U_1 = 1 #initial condition
no_simulations = 10000 # specify number of simulations wanted
target_skew = 1
possible_resistors = []
a = find_k(target_skew) #skewness param
for i in range(10000):
    possible_resistors = np.append(possible_resistors,np.random.gamma(a))
#print("possible res",possible_resistors)
actual_skew = stats.skew(possible_resistors)
print("actual skew",actual_skew)

def resistor_distribution():
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
    V = np.matrix(V)    
    H = np.matrix(H)
    ##############################################################################
    return V, H

def simulate():
    V, H = random_generation()
    A = V
    i = 1
    invL = np.array([])
    sigma = np.array([])
    #all values of 1/L and sigma_N are stored in those arrays for plotting
    while i < iteration:
        i += 1
        V, H = random_generation()
        A = V + A * np.linalg.inv(I+ H*A)
        invA = np.linalg.inv(A)
        sigma_N = 1/invA[0,0]
        sigma_N = sigma_N / i
        invL = np.append(invL, [1/i])
        sigma = np.append(sigma, [sigma_N])
    
    #find resistance of network
    effective_R = sigma_N * N 
    effective_R = 1 / effective_R
    #for justification on this calc for effective resistance
    #see note book p14
    return effective_R

x_val=np.array([])
y_val=np.array([])
z_val=np.array([])

i = 0
fig = plt.figure()
ax = plt.axes(projection='3d')
plt.xlabel("skewness")
plt.ylabel("N")
while i < 10:
    while target_skew < 5:
        simulation_results = []
        j = 1
        while j <= no_simulations:
            effective_R = simulate()
            simulation_results = np.append(simulation_results, [effective_R])
            j = j+1
            
        #find expected effective resistance
        expected_R = np.mean(possible_resistors)
        expected_R = np.exp(expected_R)
        print("exp R",expected_R)
        #dykhne resistance is just expected value of the distribution of log of 
        #resistors, see page 15 for proof
        average_R = np.mean(simulation_results)
        print("avg eff r",average_R)
        std = np.std(simulation_results)
        diff = ((expected_R - stats.describe(simulation_results).mean)/stats.describe(simulation_results).mean)*100
        diff = diff
        print("% err",diff)
        #percentage deviation
        x_val = np.append(x_val, [actual_skew])
        y_val = np.append(y_val, [diff])
        z_val = np.append(z_val, [N])
        target_skew = target_skew+0.5
        print("new target skew", target_skew)
        if target_skew < 5:
            possible_resistors = []
            a = find_k(target_skew) #skewness param
            for j in range(10000):
                possible_resistors = np.append(possible_resistors,np.random.gamma(a))
            actual_skew = stats.skew(possible_resistors)
            print("actual skew",actual_skew)

    ax.plot3D(x_val,z_val,y_val, 'red')
    target_skew = 1
    print("target skew",target_skew)
    possible_resistors = []
    a = find_k(target_skew) #skewness param
    for j in range(10000):
        possible_resistors = np.append(possible_resistors,np.random.gamma(a))
    actual_skew = stats.skew(possible_resistors)
    print("actual skew",actual_skew)
    N = N+5
    iteration = N
    I = np.identity(N)
    x_val=np.array([])
    y_val=np.array([])
    z_val=np.array([])
    i = i + 1
