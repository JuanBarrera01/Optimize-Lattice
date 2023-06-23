from math import inf
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

adj=sp.sparse.lil_array((n**2,n**2), dtype=bool)

#list of infected sites
infection_list=np.array(2*n*n, dtype=int)

#Creating a p=0 lattice with periodic boundaries
for j in range(n**2):
    randomizer = np.random.rand(3)
    #bonds with right neighbour
    #randomizer removes bond with probability p
    if j % n != (n-1) and randomizer[0]>p:
        adj[j, j+1]=adj[j+1, j]=1
    #bonds horizontal periodic
    elif randomizer[0]>p:
        adj[j, j-(n-1)]=adj[j-(n-1), j]=1

    if j < (n*(n-1)) and randomizer[1]>p:
        adj[j, j+n]=adj[j+n,j]=1
    elif randomizer[1]>p:
        adj[j,j-n*(n-1)]=adj[j-n*(n-1),j]=1
    #infection?
    #if randomizer[2] < lam/(1+lam):
    infection_list=np.append(infection_list, j)  
    
#print (infection_list)
pinfected = np.array(infection_list.size-1, dtype=int)


#looping many timesteps
for k in range (timesteps):
    #one timestep in every one of these if statements
    if (infection_list.size > 1):
        #run n_infected Monte Carlo Steps
        for j in prange(infection_list.size-1):
            #choose random infected node
            if (infection_list.size !=2):
                i=np.random.randint(1, infection_list.size-1)
            else:
                i=1
            randomizer= np.random.rand(1)
            #determining infection or heal
            if  1/(1+lam)>randomizer:
                infection_list[i]=infection_list[-1]
                infection_list=infection_list[:-1]
                #print ('heal')

            else: 
                #here we list the neighbors
                #print ('Inspected Infected Site ' + str(i))
                c=adj.getcol(i)
                neighcount = c.getnnz()
                c=c.nonzero()
                c=c[0]
                #here we choose the neighbor
                if c.size!= 0:
                    infection = np.random.choice(c)
                    if infection not in infection_list:
                        infection_list=np.append(infection_list, infection)
                        #print ('infect site ' + str(infection))
                    #else: print('already infected ' +str(infection) )
                #else: print('no friends')

            #print ( 'Cycle ' + str(k) + ' Infected after ' + str(j))
            #print (infection_list)
    pinfected = np.append(pinfected, infection_list.size-1)
    print(pinfected)
    #print (str(k) + 'steps')

#graphing pinfected
fig = plt.figure(figsize = (5,5))
plt.plot(pinfected, c = "red", marker = ".", linestyle = "None", label = "Infected", markersize = 3)
plt.xlabel("t (MC steps)", size = 12)
plt.ylabel("counts", size = 12)
plt.legend(fontsize = 12, markerscale = 2)
plt.grid(alpha = 0.5, linestyle = "--")
fig.savefig("one run.pdf", transparent = True, bbox_inches = 'tight', pad_inches = 0)

#The next function would help me count the number of bonds left. My method produces good enough results for large n
#print (a.getnnz())

#convert to nx
#G = nx.from_scipy_sparse_array(a)


#draw the graph
#nx.draw(G)
#plt.savefig("path.png")
