import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.stats import skewnorm
from mpl_toolkits import mplot3d
from scipy.interpolate import interp1d
from scipy.stats import gamma
from scipy.stats import lognorm
import random

state = 1

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

N = 25 #define the dimension of the matrix
iteration = N #no. of iterations or columns for each simulation
I = np.identity(N)
U_1 = 1 #initial condition
no_simulations = 1000 # specify number of simulations wanted
target_skew = -0.9
a = a_corresp(target_skew) #skewness param
nstep = 1000
xmin = 0
xmax = np.pi/4
def f(x):
    if state == 4:
        return np.tan(x)
    elif state == 5:
        return abs(1/(1-x))
    elif state == 6:
        return x*(x-1)*(x-2)
    elif state == 7:
        return abs((x**4)*np.log(x))
    elif state == 8:
        return 2.5*np.sin(x**2)
    elif state == 9:
        return (x - 1) * (x + 3)**2 * (x - 2)**3 * (x + 1)**4
    else:
        return 1
    
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
    #skewnormal distribution
    if state == 1:
        resistor = skewnorm.rvs(a, size=1)
    elif state == 2:
        resistor = gamma.rvs(a, scale = 1.0/a, size=1)
    elif state == 3:
        resistor = lognorm.rvs(a, size=1)
    #values here are taken as log of the resistor value
    else:
        resistor = random.choice(possible_resistors)
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

skew_plt = []
err_plt = []
while target_skew <= 0.9:
    simulation_results = []
    j = 1
    while j <= no_simulations:
        effective_R = simulate()
        simulation_results = np.append(simulation_results, [effective_R])
        j = j+1
            
    #find expected effective resistance
    expected_R = np.exp(skewnorm.expect(args=(a,)))
    #nobs = no of observations
    diff = ((expected_R - stats.describe(simulation_results).mean)/stats.describe(simulation_results).mean)*100
    #percentage deviation
    skew_plt = np.append(skew_plt, [target_skew])
    err_plt = np.append(err_plt, [diff])
    target_skew = target_skew+0.1
    if target_skew <= 0.9:
        a = a_corresp(target_skew)

state = 2
def find_k(skewness):
    k = 4.0 / (skewness**2)
    return k
target_skew = 0.1
a = find_k(target_skew)
skew_plt2 = []
err_plt2 = []
while target_skew <= 0.9:
    simulation_results = []
    j = 1
    while j <= no_simulations:
        effective_R = simulate()
        simulation_results = np.append(simulation_results, [effective_R])
        j = j+1
            
    expected_R = np.exp(gamma.mean(a,scale=1.0/a))
    diff = ((expected_R - stats.describe(simulation_results).mean)/stats.describe(simulation_results).mean)*100
    #percentage deviation
    skew_plt2 = np.append(skew_plt2, [target_skew])
    err_plt2 = np.append(err_plt2, [diff])
    target_skew = target_skew+0.1
    if target_skew <= 0.9:
        a = find_k(target_skew)

state = 3
#####################################################################
skewness_lognorm = []
a_param_lognorm = []
a = 0.0001
while a <= 0.5:
    mean, var, skew, kurt = lognorm.stats(a, moments='mvsk')
    skewness_lognorm = np.append(skewness_lognorm, skew)
    a_param_lognorm= np.append(a_param_lognorm, a)
    a = a+0.0001
b_corresp = interp1d(skewness_lognorm,a_param_lognorm)
#linear interpolation. given a skew, find corresponding a_param
#####################################################################

#plt.plot(a_param_lognorm,skewness_lognorm)
target_skew = 0.1
a = b_corresp(target_skew)
skew_plt3 = []
err_plt3 = []
while target_skew <= 0.9:
    simulation_results = []
    j = 1
    while j <= no_simulations:
        effective_R = simulate()
        simulation_results = np.append(simulation_results, [effective_R])
        j = j+1
            
    #find expected effective resistance
    expected_R = np.exp(lognorm.expect(args=(a,)))
    #nobs = no of observations
    diff = ((expected_R - stats.describe(simulation_results).mean)/stats.describe(simulation_results).mean)*100
    #percentage deviation
    skew_plt3 = np.append(skew_plt3, [target_skew])
    err_plt3 = np.append(err_plt3, [diff])
    target_skew = target_skew+0.1
    if target_skew <= 0.9:
        a = b_corresp(target_skew)

state = 4
xmin = 0
xmax = np.pi/4
possible_resistors = prob_list_gen(xmin,xmax,nstep)
list_of_percentage_error = []
for i in range(1000):
    effective_R = simulate()
    predicted_R = sum(possible_resistors) / len(possible_resistors)
    predicted_R = np.exp(predicted_R)
    percentage_deviation = ((effective_R-predicted_R)/predicted_R)*100
    list_of_percentage_error = np.append(list_of_percentage_error,percentage_deviation)
avg_error4 = sum(list_of_percentage_error)/len(list_of_percentage_error)
skew4 = stats.skew(possible_resistors)

state = 5
xmin = 0
xmax = 0.8
possible_resistors = prob_list_gen(xmin,xmax,nstep)
list_of_percentage_error = []
for i in range(1000):
    effective_R = simulate()
    predicted_R = sum(possible_resistors) / len(possible_resistors)
    predicted_R = np.exp(predicted_R)
    percentage_deviation = ((effective_R-predicted_R)/predicted_R)*100
    list_of_percentage_error = np.append(list_of_percentage_error,percentage_deviation)
avg_error5 = sum(list_of_percentage_error)/len(list_of_percentage_error)
skew5 = stats.skew(possible_resistors)

state = 6
xmin = 2
xmax = 3
possible_resistors = prob_list_gen(xmin,xmax,nstep)
list_of_percentage_error = []
for i in range(1000):
    effective_R = simulate()
    predicted_R = sum(possible_resistors) / len(possible_resistors)
    predicted_R = np.exp(predicted_R)
    percentage_deviation = ((effective_R-predicted_R)/predicted_R)*100
    list_of_percentage_error = np.append(list_of_percentage_error,percentage_deviation)
avg_error6 = sum(list_of_percentage_error)/len(list_of_percentage_error)
skew6 = stats.skew(possible_resistors)

state = 7
xmin = 0.2
xmax = 0.8
possible_resistors = prob_list_gen(xmin,xmax,nstep)
list_of_percentage_error = []
for i in range(1000):
    effective_R = simulate()
    predicted_R = sum(possible_resistors) / len(possible_resistors)
    predicted_R = np.exp(predicted_R)
    percentage_deviation = ((effective_R-predicted_R)/predicted_R)*100
    list_of_percentage_error = np.append(list_of_percentage_error,percentage_deviation)
avg_error7 = sum(list_of_percentage_error)/len(list_of_percentage_error)
skew7 = stats.skew(possible_resistors)

state = 8
xmin = 0
xmax = 1.8
possible_resistors = prob_list_gen(xmin,xmax,nstep)
list_of_percentage_error = []
for i in range(1000):
    effective_R = simulate()
    predicted_R = sum(possible_resistors) / len(possible_resistors)
    predicted_R = np.exp(predicted_R)
    percentage_deviation = ((effective_R-predicted_R)/predicted_R)*100
    list_of_percentage_error = np.append(list_of_percentage_error,percentage_deviation)
avg_error8 = sum(list_of_percentage_error)/len(list_of_percentage_error)
skew8 = stats.skew(possible_resistors)

state = 9
xmin = -0.5
xmax = 1
possible_resistors = prob_list_gen(xmin,xmax,nstep)
list_of_percentage_error = []
for i in range(1000):
    effective_R = simulate()
    predicted_R = sum(possible_resistors) / len(possible_resistors)
    predicted_R = np.exp(predicted_R)
    percentage_deviation = ((effective_R-predicted_R)/predicted_R)*100
    list_of_percentage_error = np.append(list_of_percentage_error,percentage_deviation)
avg_error9 = sum(list_of_percentage_error)/len(list_of_percentage_error)
skew9 = stats.skew(possible_resistors)

plt.plot(skew_plt,err_plt,label='skewnorm')
plt.plot(skew_plt2,err_plt2,label='gamma')
plt.plot(skew_plt3,err_plt3,label='lognorm')
plt.plot(skew4,avg_error4,'x',label='tan(x),0<x<pi/4')
plt.plot(skew5,avg_error5,'x',label='1/(1-x),0<x<0.8')
plt.plot(skew6,avg_error6,'x',label='x*(x-1)*(x-2),2<x<3')
plt.plot(skew7,avg_error7,'x',label='|(x^4)*log(x)|,0.2<x<0.8')
plt.plot(skew8,avg_error8,'x',label='2.5*sin(x^2),0<x<1.8')
plt.plot(skew9,avg_error9,'x',label='(x-1)*(x+3)^2*(x-2)^3*(x+1)^4,-0.5<x<1')
plt.xlabel("skewness")
plt.legend(loc='upper center', bbox_to_anchor=(0.5, -0.05), shadow=True, ncol=2)
plt.ylabel("% error")
plt.title("25 by 25 networks, unimodal distributions with |skewness| < 1")
