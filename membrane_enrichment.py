#!/usr/bin/env python3
#
#  membrane_enrichment.py: Calculates local enrichment of for each lipid 
#               species in each leaflet of a flat bilayer membrane. What is 
#               local is defined by bin size.
#               Note: The program assumes that the membrane is composed of 2 
#               lipid species and cholesterol.
#  Akshara Sharma

# Sample command to run:
# python_call program_name system_file trajectory_file sel_string_1 sel_string_2 sel_string_3 grid_size histogram_bins 
# python3 membrane_enrichment.py 1DIPC.psf 1DIPC_8us.dcd 'resname=="DIPC"' 'resname=="DPPC"' 'resname=="CHOL"' 10 10

import os
import sys
import loos
import loos.pyloos
import loosfunc
import math
import numpy as np


dir = os.path.dirname(__file__)

#header = " ".join(sys.argv)
#print("# ", header)

system_file = sys.argv[1]
traj_file = sys.argv[2]
lipid1_sel_string = sys.argv[3]
lipid2_sel_string = sys.argv[4]
chol_sel_string = sys.argv[5]
num_bins = int(sys.argv[6])
hist_bins = int(sys.argv[7])

system = loos.createSystem(system_file)
traj = loos.pyloos.Trajectory(traj_file, system)
sel1 = loos.selectAtoms(system, lipid1_sel_string)
sel2 = loos.selectAtoms(system, lipid2_sel_string)
selc = loos.selectAtoms(system, chol_sel_string)

# break sel1, sel2, and selc into individual lipid molecule selections
selections = [sel1.splitByMolecule(),sel2.splitByMolecule(),selc.splitByMolecule()]

# create empty arrays to which histogram will be appended for each species in each leaflet.
hist = [[np.zeros(shape=[0, 2]) for i in range(len(selections))] for j in range(2)]

frame = 0
while traj.nextFrame():
    # break species selections into upper, lower leaflets
    # index 0 is upper, index 1 is lower leaflet
    sel1_upper, sel1_lower = loosfunc.leafletLipidSeparator(selections[0])
    sel2_upper, sel2_lower = loosfunc.leafletLipidSeparator(selections[1])
    selc_upper, selc_lower = loosfunc.leafletLipidSeparator(selections[2])
    leaflet_sel = [[sel1_upper,sel2_upper,selc_upper],[sel1_lower,sel2_lower,selc_lower]]
    
    mol_count = np.array([[len(leaflet_sel[0][0]), len(leaflet_sel[0][1]), len(leaflet_sel[0][2])],
                          [len(leaflet_sel[1][0]), len(leaflet_sel[1][1]), len(leaflet_sel[1][2])]])
    
    size = system.periodicBox()
    bin_size = size/num_bins
    
    # iterate over leaflets. index 0 is upper, index 1 is lower leaflet
    for l in range(2):
        # calculate 2D ideal gas concentration for each species in leaflet
        ideal_conc = mol_count[l]/(size[0]*size[1])
        
        # calculate number of sub-bin shifts and shift distance
        n_incr = math.ceil( np.sqrt(np.sum(mol_count[l]))/num_bins )
        incr = bin_size/n_incr
        
        # loop over species in leaflet
        for i in range(len(leaflet_sel[l])):
            enr = np.zeros((0,num_bins))
            
            # iterating over useful range of grid shifts
            for y in range(n_incr):
                for x in range(n_incr):
                    shift = np.array((incr[0], 0.0, 0.0))
                    if x == 0:
                        shift = shift + np.array((0.0, incr[1], 0.0))
                    Gshift = loos.GCoord(*shift)
                    
                    # loop over lipids, assign to leaflet, bin
                    pop = np.zeros((num_bins,num_bins))
                    for element in leaflet_sel[l][i]:
                        element.translate(Gshift)
                        element.reimageByAtom()
                        for atom in element:
                            loosfunc.assign_bin_atom(atom, bin_size, pop)
                    
                    # calculate fractional number of lipids in bin
                    pop /= len(leaflet_sel[l][i][0])
                    
                    # iterate over bins, calculate enrichment, append to enrichment array
                    enr = np.vstack((enr, loosfunc.enrichment(pop, bin_size, ideal_conc[i])))
            
            # calculate and append histogams of enrichment to the histogram column stack
            labels = [sel1[0].resname(),sel2[0].resname(),selc[0].resname()]
            leaflet_label = ["Upper","Lower"]
            np.savetxt(os.path.join(dir,"CV Data",leaflet_label[l]+'_leaflet_enrichment_'+labels[i]+'_'+str(frame)+'.txt'), enr.flatten(), fmt='%.6e')
            hist[l][i] = np.vstack((hist[l][i], loosfunc.histogram(enr, hist_bins, 0.0, 3.0)))
    frame += 1

# output histograms of each lipid species in each leaflet
labels = [sel1[0].resname(),sel2[0].resname(),selc[0].resname()]
leaflet_label = ["Upper","Lower"]
for l in range(2):
    for i in range(len(leaflet_sel[l])):
        np.savetxt(os.path.join(dir,"CV Data",leaflet_label[l]+'_leaflet_enrichment_'+labels[i]+'.txt'), hist[l][i], fmt='%.6e')
