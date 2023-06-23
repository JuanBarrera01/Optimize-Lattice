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
print('nruns')
nruns=input()
nruns=int(nruns)
import matplotlib.colors as mcolors



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
        trace = np.array(1, dtype=int)
        t=0
    return lattice, inflist, trace, t

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
                    elif infection==1:
                        if inflist[i] < (size*(size-1)):
                           if inflist[i]+size not in inflist:
                                inflist=np.append(inflist, inflist[i]+size)
                        elif inflist[i]-size*(size-1) not in inflist:
                            inflist=np.append(inflist, inflist[i]-size*(size-1))
                    elif infection==2:
                        if inflist[i]% size != 0: 
                            if inflist[i]-1 not in inflist:
                                inflist=np.append(inflist, inflist[i]-1)
                        elif inflist[i]+size-1 not in inflist:
                            inflist=np.append(inflist, inflist[i]+size-1)
                    elif infection==3:
                        if inflist[i]>size-1: 
                            if inflist[i]-size not in inflist:
                                inflist=np.append(inflist, inflist[i]-size)
                        elif inflist[i]+size*(size-1) not in inflist:
                            inflist=np.append(inflist, inflist[i]+size)
                    #not in inflist:
                    #inflist=np.append(inflist, infection)
                        #print ('infect site ' + str(infection))
                    #else: print('already infected ' +str(infection) )
            t+=1/inflist.size
            if t >= trace.size :
                trace = np.append(trace,(inflist.size-1)/(size**2))
    else:
        trace = np.append(trace,(inflist.size-1)/(size**2))
    return inflist, trace, t

def cycle(inflist, infrate, lattice, size, trace, t, duration):
    for k in range(duration):
        inflist, trace, t=step(inflist, infrate, lattice, size , trace, t)
        print(k)
    return  trace

def simulation(size, perc, infrate, duration, numruns):
    s = {}
    tlength = np.zeros(numruns)
    for trial in prange(numruns):
        lattice, inflist, trace, t = create_lattice(size, perc)
        trace = cycle(inflist, infrate, lattice, size, trace, t, duration)
        s[trial]=trace
        tlength[trial]=trace.size
    return s, tlength

def plotmany(manyruns, numruns):
    fig = plt.figure(figsize = (5,5))
    colors=list(mcolors.TABLEAU_COLORS)
    for trial in range (numruns):
        plt.plot( manyruns[trial] , c = colors[trial], marker = ".", linestyle = "None", label = "trial " + str(trial), markersize = 3)
    plt.xlabel("t (MC steps)", size = 12)
    plt.ylabel("counts", size = 12)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(fontsize = 12, markerscale = 2)
    plt.grid(alpha = 0.5, linestyle = "--")
    plt.title("Size " + str(n) + "p= " + str(p) + "lambda " + str(lam))
    fig.savefig("Graphs/Multirun/" "Size " + str(n) + "p= " + str(p) + "lambda " + str(lam) + "steps " + str(timesteps) + str(numruns) + " runs.jpg", transparent = True, bbox_inches = 'tight', pad_inches = 0)
    fig.savefig("Graphs/Multirun/" "Size " + str(n) + "p= " + str(p) + "lambda " + str(lam) + "steps " + str(timesteps) + str(numruns) + " runs.pdf", transparent = True, bbox_inches = 'tight', pad_inches = 0)

def plotav(manyruns, size, infrate, numruns):
    shortest=min(size)
    shortest=int(shortest)
    average = np.zeros(shortest)
    for i in range (shortest):
        for j in range(numruns):
            average[i]+=manyruns[j][i]
    average=average/numruns
    fig = plt.figure(figsize = (5,5))
    plt.plot( average , c = 'red' , marker = ".", linestyle = "None", label = infrate, markersize = 3)
    plt.xlabel("t (MC steps)", size = 12)
    plt.ylabel("counts", size = 12)
    plt.xscale("log")
    plt.yscale("log")
    plt.legend(fontsize = 12, markerscale = 2)
    plt.grid(alpha = 0.5, linestyle = "--")
    plt.title("Size " + str(n) + "p= " + str(p) + "lambda " + str(lam))
    fig.savefig(  "Graphs/Averages/" "Size " + str(n) + "p= " + str(p) + "lambda " + str(lam) + "steps " + str(timesteps) + " average " + str(numruns) + " runs.jpg", transparent = True, bbox_inches = 'tight', pad_inches = 0)
    fig.savefig(  "Graphs/Averages/" "Size " + str(n) + "p= " + str(p) + "lambda " + str(lam) + "steps " + str(timesteps) + " average " + str(numruns) + " runs.pdf", transparent = True, bbox_inches = 'tight', pad_inches = 0)

samples, triallength = simulation(n, p, lam, timesteps, nruns)

plotmany(samples, nruns)

plotav(samples, triallength, lam, nruns)

#The next function would help me count the number of bonds left. My method produces good enough results for large n
#print (a.getnnz())

#convert to nx
#G = nx.from_scipy_sparse_array(a)


#draw the graph
#nx.draw(G)
#plt.savefig("path.png")






