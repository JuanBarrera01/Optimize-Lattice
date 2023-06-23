from math import inf
from os import times
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import scipy as sp
from numpy import linalg
from numpy.random import rand
from numba import jit, prange

#autopilot
#size
print ('size')
n=input()
n = int(n)
#percolation
print('percolation')
p = input()
p= float(p)
print('lambda')
lam = input()
lam = float (lam)
print('timesteps')
timesteps = input()
timesteps = int(timesteps)




#list of infected sites


#Creating a lattice with periodic boundaries
def create_lattice(size, perc):
    lattice=np.zeros((size**2,4), dtype=bool)
    #use order right, down, left, up
    inflist=np.array(2*size**2, dtype=int)
    for j in range(size**2):
        randomizer = np.random.rand(3)
        #bonds with right neighbour
        #randomizer removes bond with probability p
        if j % size != (size-1) and randomizer[0]>=perc:
            lattice[j, 0]=lattice[j+1, 2]=1
        #bonds horizontal periodic
        elif randomizer[0]>=perc:
            lattice[j, 0]=lattice[j-(size-1), 2]=1
        if j < (size*(size-1)) and randomizer[1]>=perc:
            lattice[j, 1]=lattice[j+size,3]=1
        elif randomizer[1]>=perc:
            lattice[j,1]=lattice[j-size*(size-1),3]=1
        #infection?
        #if randomizer[2] < lam/(1+lam):
        inflist=np.append(inflist, j)  
        trace = np.array([0, 1], dtype=float)
        t=0
    return lattice, inflist, trace, t


adj, infection_list, ninfected, time = create_lattice(n, p)
print('matrix ready')


def step(inflist, infrate, lattice, size, trace, t):
    if (inflist.size > 1):
        #run n_infected Monte Carlo Steps
        for j in prange(inflist.size-1):
            #choose random infected node
            if (inflist.size !=2):
                i=np.random.randint(1, inflist.size-1)
            else:
                i=1
            randomizer= np.random.rand(1)
            #determining infection or heal
            if  1/(1+infrate)>randomizer:
                inflist[i]=inflist[-1]
                inflist=inflist[:-1]

            else:
                infection = np.random.randint(0,3)
                if lattice[inflist[i]][infection]== 1:
                    if infection==0:
                        if inflist[i] % size != (size-1):
                           if inflist[i]+1 not in inflist:
                               inflist=np.append(inflist, inflist[i]+1)
                        elif inflist[i]-(size-1) not in inflist:
                            inflist=np.append(inflist, inflist[i]-(size-1))
                    if infection==1:
                        if inflist[i] < (size*(size-1)):
                           if inflist[i]+size not in inflist:
                                inflist=np.append(inflist, inflist[i]+size)
                        elif inflist[i]-size*(size-1) not in inflist:
                            inflist=np.append(inflist, inflist[i]-size*(size-1))
                    if infection==2:
                        if inflist[i]% size != 0: 
                            if inflist[i]-1 not in inflist:
                                inflist=np.append(inflist, inflist[i]-1)
                        elif inflist[i]+size-1 not in inflist:
                            inflist=np.append(inflist, inflist[i]+size-1)
                    if infection==3:
                        if inflist[i]>size-1: 
                            if inflist[i]-size not in inflist:
                                inflist=np.append(inflist, inflist[i]-size)
                        elif inflist[i]+size*(size-1) not in inflist:
                            inflist=np.append(inflist, inflist[i]+size)
                    #not in inflist:
                    #inflist=np.append(inflist, infection)
                        #print ('infect site ' + str(infection))
                    #else: print('already infected ' +str(infection) )
            trace = np.append(trace, [float(trace[2*t]) + 1/inflist.size , (inflist.size-1)/(size**2)])
            t+=1
    return inflist, trace, t

for k in range(timesteps):
    print(k)
    infection_list, ninfected, time=step(infection_list, lam, adj, n , ninfected, time)


print(infection_list)


ninfected=np.reshape(ninfected, ( int(ninfected.size/2), 2))

fig = plt.figure(figsize = (5,5))
plt.plot(ninfected[:, 0], ninfected[:, 1] , c = "red", marker = ".", linestyle = "None", label = "Infected", markersize = 3)
plt.xlabel("t (MC steps)", size = 12)
plt.ylabel("counts", size = 12)
plt.legend(fontsize = 12, markerscale = 2)
plt.grid(alpha = 0.5, linestyle = "--")
fig.savefig("Size " + str(n) + "p= " + str(p) + "lambda " + str(lam) + "steps " + str(timesteps) + " one run.jpg", transparent = True, bbox_inches = 'tight', pad_inches = 0)

#The next function would help me count the number of bonds left. My method produces good enough results for large n
#print (a.getnnz())

#convert to nx
#G = nx.from_scipy_sparse_array(a)


#draw the graph
#nx.draw(G)
#plt.savefig("path.png")






