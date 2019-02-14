import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import skewnorm
from mpl_toolkits import mplot3d
from scipy.interpolate import interp1d

#this code is to test dykhne thm using skewnorm distribution

#####################################################################
skewness = []
a_param = []
a = -10
while a <= 10:
    mean, var, skew, kurt = skewnorm.stats(a, moments='mvsk')
    skewness = np.append(skewness, skew)
    a_param= np.append(a_param, a)
    a = a+0.001
a_corresp = interp1d(skewness,a_param)
#linear interpolation. given a skew, find corresponding a_param
#####################################################################

N = 5 #define the dimension of the matrix
iteration = N #no. of iterations or columns for each simulation
I = np.identity(N)
U_1 = 1 #initial condition
no_simulations = 1000 # specify number of simulations wanted
target_skew = -0.95
a = a_corresp(target_skew) #skewness param

def multiplyList(numbers):  
    #exponentiate elements in a list and multiply them together
    #e.g. numbers = (ln(1),ln(3),ln(5)) gives 15
    total = 1
    for x in numbers:
        total *= np.exp(x)  
    return total  

def resistor_distribution():
    #define your random distribution here
    #skewnormal distribution
    resistor = skewnorm.rvs(a, size=1)
    #values here are taken as log of the resistor value
    resistor = np.exp(resistor)
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
    while target_skew < 1:
        simulation_results = []
        j = 1
        while j <= no_simulations:
            effective_R = simulate()
            simulation_results = np.append(simulation_results, [effective_R])
            j = j+1
            
        #find expected effective resistance
        expected_R = np.exp(skewnorm.expect(args=(a,)))
        #dykhne resistance is just expected value of the distribution of log of 
        #resistors, see page 15 for proof
        average_R = np.mean(simulation_results)
        std = np.std(simulation_results)

        #nobs = no of observations
        diff = ((expected_R - stats.describe(simulation_results).mean)/stats.describe(simulation_results).mean)*100
        diff = abs(diff)
        #percentage deviation
        x_val = np.append(x_val, [target_skew])
        y_val = np.append(y_val, [diff])
        z_val = np.append(z_val, [N])
        target_skew = target_skew+0.1
        if target_skew < 1:
            a = a_corresp(target_skew)

    ax.plot3D(x_val,z_val,y_val, 'red')
    target_skew = -0.95
    a = a_corresp(target_skew)
    N = N+5
    iteration = N
    I = np.identity(N)
    x_val=np.array([])
    y_val=np.array([])
    z_val=np.array([])
    i = i + 1
