import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import random

'''
things to vary:
line 11, N paramter, size of network
line 20, distribution of network (currently set to uniform between 1 to 10)
graphs are outputed to same location as where the code is in jpg format
'''
N = 4 #define the dimension of the matrix
original_N = N
iteration = N #no. of iterations or columns for each simulation
I = np.identity(N)
U_1 = 1

def resistor_distribution():    
    #define your random distribution here
    #resistor = np.exp(random.choice(possible_resistors))
    resistor = np.random.randint(1, 10)
    #resistor = 1
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
        #H[i,i] = 999.0
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
    return V_resistors, H_resistors, A, V_list, H_list, effective_R

V_resistors, H_resistors, A, V_list, H_list, effective_R = simulate()
#find expected effective resistance
#print("effective R is: ", effective_R)

#solve simultaneous equations to find I_1 & U_i
#-A_11U_1 = -I_1 + A_12U_2 + A_13U_3
#-A_21U_1 = 0I_1 + A_22U_2 + A_23U_3
#...
print("eff r is ", effective_R)
pred_R = 5.0
perc_err = abs(pred_R-effective_R)/effective_R*100
print("percentage error", perc_err)

#print("now we try to find I_1 & U_i \n")
A = np.array(A)
a = np.array(A) #matrix of coeff of unknowns
b = np.array([]) #sol of simultaneous eqn at lhs
i = 0
U = [U_1]
I = np.zeros(N)

while i < N :
    a[i][0] = 0
    b = np.append(b, [-U_1 * A[i][0]])
    i+=1
a[0][0] = -1
sol = np.linalg.solve(a,b)
#print("I_1 is \n", sol[0])
I[0] = sol[0]
i = 1
while i < N :
    #print("U_", i+1 ,"is \n", sol[i])
    U = np.append(U,sol[i])
    i += 1

#We now have U, I, A, V_list, H_list
#print("A_N is: ", A)
print("I_N is: ", I)
print("U_N is: ", U)
U = np.array(U)
I = np.array(I)
U_list = [U.reshape((N,1))]*N
I_horizontal_list = [np.ones((N,1)) for i in range(N)]
I_vertical_list = [np.ones((N,1)) for i in range(N)]

print("V_list\n",V_list)
print("H_list\n",H_list)
N=original_N
N=N-1

while N >= 0:
    if N==original_N-1:
        I_horizontal_list[N]=I.reshape((original_N,1))
        #print("first time. I_horizontal_",N,"is: ",I_horizontal_list[N])
    else:
        #I_horizontal_list[N] = np.array(I_horizontal_list[N+1] - V_list[N+1]*U_list[N+1])
        I_horizontal_list[N] = I_horizontal_list[N+1] - np.matrix(V_list[N+1])*np.matrix(U_list[N+1])
        #print("A is",V_list[N+1])
        #print("B is",U_list[N+1])
        #print("A*B is",np.dot(V_list[N+1],U_list[N+1]))
        '''
        for m in range(original_N):
            for n in range(original_N):
                if I_horizontal_list[m][n] < 10**(-10):
                    I_horizontal_list[m][n] = 0.0
                elif I_horizontal_list[m][n] > I_horizontal_list[original_N-1][0]:
                    I_horizontal_list[m][n] = I_horizontal_list[original_N-1][0]
        '''
        #print("I_horizontal_",N,"is: ",I_horizontal_list[N])
    if N!=original_N-1:        
        U_list[N] = U_list[N+1] - H_list[N+1]*I_horizontal_list[N] 
        '''
        for m in range(original_N):
            for n in range(original_N):
                if U_list[m][n] < 10**(-10):
                    U_list[m][n] = 0.0
        '''
        #print("U_list",N,"is: \n",U_list[N])
    j = 0
    while j < (original_N-1):
        I_vertical_list[N][j] = (U_list[N][j]-U_list[N][j+1])/V_resistors[N][j]
        j = j + 1
    I_vertical_list[N][j] = U_list[N][j]/V_resistors[N][j]
    '''
    for m in range(original_N):
        for n in range(original_N):
            if I_vertical_list[m][n] > I_horizontal_list[original_N-1][0]:
                I_vertical_list[m][n] = I_horizontal_list[original_N-1][0]
    '''
    #print("I_vertical_",N,"is: ",I_vertical_list[N])
    #input("Press Enter to continue...")
    N = N-1
    
i = 0
sum_vert = 0
while i < original_N:
    sum_vert = sum_vert + I_vertical_list[i][0]
    i = i+1
print("sum of first row of vertical current is",sum_vert)
print("top right entry is",I_horizontal_list[original_N-1][0])
deviation_in_current = I_horizontal_list[original_N-1][0]-sum_vert
print("deviation",deviation_in_current)

potential_drop = []
for m in range(original_N):
    for n in range(original_N):
        drops = I_vertical_list[m][n]*V_resistors[m][n]
        potential_drop = np.append(potential_drop,drops)
pd_error = sum(potential_drop)/original_N
pd_error = abs(pd_error-U_1)/pd_error*100
print("average % error on potential drops across each column is: ",pd_error)

V_resistors = np.round(V_resistors,2)
H_resistors = np.round(H_resistors,2)
#print("vertical resistors:",V_resistors)
#print("horizontal resistors:",H_resistors)

plt.figure()
N=original_N
G=nx.DiGraph()
i=1
j=N+1
k=1
while i<=N+1:
    while j>0:
        G.add_node(k, pos=(i,j))
        k+=1
        j-=1
    j=N+1
    i+=1 
G.remove_node(k-1)

i=1
while i<=N:
    j=0
    while j<N:
        G.add_edge(i+j*(N+1),i+(j+1)*(N+1),weight=H_resistors[j][i-1])
        j+=1
    i+=1
    
i=1
while i<=N:
    j=0
    while j<N:
        G.add_edge(i+j*(N+1),i+j*(N+1)+1,weight=V_resistors[j][i-1])
        j+=1
    i+=1

pos=nx.get_node_attributes(G,'pos')
nx.draw(G,pos,arrows=False)
labels = nx.get_edge_attributes(G,'weight')
nx.draw_networkx_edge_labels(G,pos,edge_labels=labels,arrows=False)
#nx.draw_networkx_nodes(G, pos)
nx.draw_networkx_labels(G, pos)
plt.title("resistor plot")
#plt.figure(figsize=(100,50))
plt.figure()

I_horizontal_list = np.array(np.round(I_horizontal_list,2))
I_vertical_list = np.array(np.round(I_vertical_list,2))
#print("horizontal currents:",I_horizontal_list)
#print("vertical currents:",I_vertical_list)

'''
#force current at row 1 to be 0 for colour contrast
j=0
while j<N:
    I_horizontal_list[j][0]=0
    j = j+1
'''

G2=nx.DiGraph()
i=1
j=N+1
k=1
while i<=N+1:
    while j>0:
        G2.add_node(k, pos=(i,j))
        k+=1
        j-=1
    j=N+1
    i+=1 
G2.remove_node(k-1)
i=1
while i<=N:
    j=0
    while j<N:
        if I_horizontal_list[j][i-1] >= 0:
            G2.add_edge(i+(j+1)*(N+1),i+j*(N+1),weight=I_horizontal_list[j][i-1])
        else:
            G2.add_edge(i+j*(N+1),i+(j+1)*(N+1),weight=I_horizontal_list[j][i-1])
        
        G2.add_edge(i+j*(N+1),i+j*(N+1)+1,weight=I_vertical_list[j][i-1])
        j+=1
    i+=1
    
pos2=nx.get_node_attributes(G2,'pos')
nx.draw(G2,pos2)
labels2 = nx.get_edge_attributes(G2,'weight')
nx.draw_networkx_edge_labels(G2,pos2,edge_labels=labels2,arrows=True)
#nx.draw_networkx_nodes(G2, pos2)
#nx.draw_networkx_labels(G2, pos2)
plt.title("current value plot")
plt.savefig("current_value.jpg")
plt.figure()

temp = list(labels2.values())
weights_list = []
for i in range(len(temp)):
    weights_list = np.append(weights_list, temp[i][0])

'''
#plot only middle 50% current
quartile1 = np.percentile(weights_list, 15)
quartile3 = np.percentile(weights_list, 85)
i = 0
while i < len(weights_list):
    if weights_list[i] < quartile1:
        weights_list[i] = quartile1
    elif weights_list[i] > quartile3:
        weights_list[i] = quartile3
    i = i+1
'''

'''
#plot only top 50% currents
for i in range(len(weights_list)):
    weights_list[i] = abs(weights_list[i])
median = np.median(weights_list)
i = 0
while i < len(weights_list):
    if weights_list[i] < median:
        weights_list[i] = 0
    i = i+1
'''

mycolors = weights_list
vmax = max(weights_list)
nodes = nx.draw_networkx_nodes(G2,pos2,node_color='k', with_labels=False,node_size=0)
edges = nx.draw_networkx_edges(G2,pos2,edge_color=mycolors,width=2, arrows=False,
                               edge_cmap=plt.cm.seismic,edge_vmax=vmax,edge_vmin=-vmax)                               
plt.colorbar(edges)
plt.axis('off')
plt.title("current plot")
plt.savefig("graph.jpg")
