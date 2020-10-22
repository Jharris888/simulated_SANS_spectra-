# JJH Debye calculation  
# updated 10.21.20 


import MDAnalysis as mda 
import MDAnalysis.analysis.rdf 
import matplotlib.pyplot as plt
import numpy as np
import math 
from itertools import combinations
from itertools import combinations_with_replacement
import MDAnalysis.analysis.distances 
import scipy.cluster
from scipy.cluster.hierarchy import single, fcluster 
import MDAnalysis.core.topologyattrs
import MDAnalysis.analysis.rdf 

################################################################
# Some important facts:                                        #
#                                                              #
# This algorithm works only for DPC micelles                   #
#                                                              #
# DPC has C, H, N, O, P atoms                                  #
#                                                              #
# f(Q) = b, the form factor of each atom is a constant         #             
#                                                              #
# I(Q) = sum(b[i]*b[j]*sin(Q*rij)/(Q*rij)), Debye              #
#                                                              #
#                                                              #
################################################################

# load the trajectories
all_u = []
for r in range(1,7):
    all_u.append(mda.Universe('reps/rep%i/tric.gro'%r,'reps/rep%i/tric.xtc'%r))
    print('loaded trajectory %i of 6'%r)


# cluster the separate regions 

bins = 5000

all_rijs = np.zeros([15,bins])
# load what we need 
for r in range(0,6):         
    u = all_u[r]
    Q = np.load('Analysis/Q_pambou.npy')

    frame = 10 # we will look at the last 'frame' frames
    bins = 5000 # how many bins to use in the histogram approximation 
    
    ####################################################################
        
        # select the clusters 
        
        # cluster things again 

    box = u.dimensions
    cutoff = 8.5 # angstroms, based on second peak of rdf graph 
    cS = u.select_atoms('name C312') # select the atoms
    r1 = [] # get positions 
    for ts in u.trajectory:
        r1.append((cS.positions))
        #print('loading positions frame %i of %i'%(u.trajectory.frame+1, len(u.trajectory)))
    r1 = np.array(r1)
    dist = [] # get all pair distances formatted as flat upper triangles
    for i in range(len(r1)):
        dist.append(MDAnalysis.analysis.distances.self_distance_array(r1[i], box=box))
        #print('loading distances frame %i of %i'%(i+1, len(u.trajectory)))
    dist = np.array(dist)
    #print('clustering.')
    z = [] # perform hierarchical single-linkage clustering 
    for i in range(len(dist)):
        z.append(single(dist[i]))
    z = np.array(z)
    #print('clustering..')
    hierarchy = [] # get clusters using cutoff (in angstroms)
    for i in range(len(z)):
        hierarchy.append(fcluster(z[i], cutoff, criterion='distance'))
    hierarchy = np.array(hierarchy) 

        # select the indices of the atoms in each cluster 
        
    DPC = u.select_atoms('resname FOS12')

    clusters = []
    for j in range(-frame,-1):
        clusters1 = []
        for i in range(1, np.amax(hierarchy[j])+1):
            c_inds = np.where(hierarchy[j] == i)
            c_sel = " ".join(list((c_inds[0]+1).astype('str')))
            sel_atoms = DPC.select_atoms('resid %s'%c_sel) # select beads in the cluster
            clusters1.append(sel_atoms)
        clusters.append(np.array(clusters1))

        # clusters == a list with an array for each frame of atom groups for each cluster 
        # find the distances rij for each type of atom
        # first we can make a list of all atom names
        
    all_heights = np.zeros([15,bins])
    for l in range(len(clusters)):
        u.trajectory[len(u.trajectory)-frame+l]
        atoms = ['name C*','name H*','name N*','name O*','name P*']
        
        box = u.dimensions
        
        # we need an array of all atom positions for each cluster 
        
        #print('finding atom positions...')
        atom_pos = []
        for i in range(len(atoms)):
            pos1 = []
            for j in range(len(clusters[l])):
                atom = clusters[l][j].select_atoms(atoms[i])
                pos1.append(atom.positions)
            atom_pos.append(np.array(pos1))
        
        
        # we can compute every possible combination of pair distances for each cluster using atom positions
        # make histogram of atom distances, we can cut out zero in this step
        # 'bins' bins will be made which have binwidth = 150/'bins'
        # binning as we go through each frame is faster and saves memory
        # something like recursively adding counts to one big histogram
    
        atom_pairs = list(combinations_with_replacement([0,1,2,3,4],2))
        
        height = []
        height1 = np.zeros(bins)
        for i in range(len(atom_pairs)):
            for j in range(len(atom_pos)):
                clust_dist = []
                for k in range(len(atom_pos[j])):
                    ref = atom_pos[atom_pairs[i][0]][k]
                    conf = atom_pos[atom_pairs[i][1]][k]
                    clust_dist1 = mda.analysis.distances.distance_array(ref,conf,box=box)
                    clust_dist.append(np.concatenate(clust_dist1))
                # concatenate the intracluster distances before making the histograms
                if len(clust_dist) > 0:
                    hist1 = np.histogram(np.concatenate(np.array(clust_dist)),bins=bins,range=(0.05,150.5))[0]
                    dists = np.histogram(np.concatenate(np.array(clust_dist)),bins=bins,range=(0.05,150.5))[1]
                    height1 = (height1 + hist1)
            height.append(height1)
        # Height is the hieght of the histogram spread over the 15 pairs of atoms
        height = np.array(height)
        all_rijs = (all_rijs+height) # all_rijs is an array for each set of pairs, for each bin: 15x5000
    print('finished rep %i'%(r+1))
np.save('Analysis/data/all_rijs.npy',all_rijs)

########################################################################################################



# define a function that can sum the intensity over a range of Q values
# rij will be computed for each pair of atoms as the center of each bin
# the sum of each rij will be multiplied by the height of each bin 
  
bins = 5000

all_rij = np.load('Analysis/data/all_rijs.npy')     

dists = np.linspace(0.05,150.5,bins) # distances used to define rij
Q = np.load('Analysis/Q_pambou.npy') # the Q scattering vector from Pambou et al. (2015)

# these are the lit. scattering lengths of each atom in DPC

Cb = 0.66
Hb = 0.65 # Hydrogen is considered as Deuterium 
Nb = 0.94
Ob = 0.58 
Pb = 0.50

# all fifteen combinations of C, H, N, O, P bij values

bij = []
bij.append(Cb*Cb)
bij.append(Cb*Hb)
bij.append(Cb*Nb)
bij.append(Cb*Ob)
bij.append(Cb*Pb)
bij.append(Hb*Hb)
bij.append(Hb*Nb)
bij.append(Hb*Ob)
bij.append(Hb*Pb)
bij.append(Nb*Nb)
bij.append(Nb*Ob)
bij.append(Nb*Pb)
bij.append(Ob*Ob)
bij.append(Ob*Pb)
bij.append(Pb*Pb)


# we are left with an array of arrays which correspond to the regions and atom pair distance heights

# function to calulate I(Q,rij)

def struc(Q,rij):
    val = (math.sin(Q*rij))/(Q*rij)
    return val

# sum over each distance [k] for each set of pairs [j] for each Q [i]  
print('summing for each Q...')
I_Q = []
for i in range(len(Q)):
    I_sum = []
    for j in range(0,15):
        I_pairs = []
        for k in range(len(dists)-1): # subtract 1 so we can find the center of each bin
            rij = (dists[k]+dists[k+1])/2 # center of each bin 
            I_pairs.append(all_rij[j][k]*bij[j]*struc(Q[i],rij)) # bij*sin(Qr)/(Qr) (Debye!)
        I_sum.append(np.sum(I_pairs))
    I_Q.append(np.sum(I_sum)) 
    print('%s of 125 sums complete..'%(i+1))

np.save('Analysis/IQs.npy',I_Q)






    