# simulated_SANS_spectra 
# Jonathan Harris
# last updated 10.21.20

 Uses the Debye scattering equation to compute the small-angle neutron scattering (SANS) profile for a surfactant system trajectory. 
 
 ################################################################
 Some important facts:                                        
                                                              
 This algorithm works only for DPC micelles                   
                                                              
 DPC has C, H, N, O, P atoms                                  
                                                              
 f(Q) = b, the form factor of each atom is a constant scattering length                  
                                                              
 I(Q) = sum(b[i]*b[j]*sin(Q*rij)/(Q*rij)), Debye scattering equation          
 
 ################################################################
                              
 1. Each frame of the input trajectories are clustered using a defined cutoff distance with single-link hierarchical clustering. The minimum image convention is applied when calculating distances. Multiple trajectories may be input at once to generate a larger ensemble average. 
 
 2. The interatomic distances (rij) are computed for the atoms within each cluster. The distances from all replicates and frames are binned into the same histograms, which are considered to be part of the same equilibrium ensemble. The pair distance histograms are separated by atom type, so that each possible pair of atom types has its own binning. The different types of pairs are separated so that the correct scattering lengths can be applied in the Debye formula. 

3. The Debye scattering function is called a separate time for each type of pair distance histogram, and for each Q scattering vector. For each type of histogram call, the correct bij scattering length constant is input to the Debye function. The output is a len(Q) list of scattering intensities, I_Q, which should be normalized by dividing by the first term, I_Q[0]

![Graphic Summary](./debye_graphic.pdf =200x200)



 
