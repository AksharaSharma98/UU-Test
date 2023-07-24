# -*- coding: utf-8 -*-
"""
Created on Mon Jan 16 16:09:08 2023

@author: saroj
"""

import numpy as np
import loos
import loos.pyloos


def categorize(selections):
    leaflet_sel = [[] for i in range(2)]
    for selection in selections:
        sel_upper, sel_lower = leafletLipidSeparator(selection)
        leaflet_sel[0].append(sel_upper)
        leaflet_sel[1].append(sel_lower)
    
    mol_count = np.zeros((2,len(selections)))
    for l in range(2):
        for s in range(len(selections)):
            mol_count[l][s] = len(leaflet_sel[l][s])
    
    return leaflet_sel, mol_count


def leafletLipidSeparator(lipids):
    '''
    Returns a tuple corresponding to lipids in Up and Lo leaflet.
    '''
    lipidUp = []
    lipidLo = []

    for lipid in lipids:
        if lipid.centroid().z() > 0:
            lipidUp.append(lipid)
        else:
            lipidLo.append(lipid)
            
    return (lipidUp, lipidLo)


def assign_bin(lipid, bin_size, pop):
    '''
    Assigns the input lipid to one bin in the n*n population array
    using the lipid centroid position.
    '''
    center = lipid.centroid()
    
    x = int(center[0]//bin_size[0])
    y = int(center[1]//bin_size[1])
    
    # assign to bin
    pop[x][y] += 1


def assign_bin_atom(atom, bin_size, pop):
    '''
    Assigns the input lipid to one bin in the n*n population array
    using the lipid centroid position.
    '''
    coords = loos.Atom.coords(atom)
    
    x = int(coords[0]//bin_size[0])
    y = int(coords[1]//bin_size[1])
    
    # assign to bin
    pop[x][y] += 1


def assign_circular_bin(lipid, shift, center, radius, pop):
    '''
    Assigns the 
    '''
    coords = lipid.centroid()

    dx = coords[0] + shift[0] - center[0]
    dy = coords[1] + shift[1] - center[1]
    
    # assign to bin
    if dx*dx + dy*dy <= radius**2:
        pop[0] += 1


def assign_circular_bin_atom(atom, shift, center, radius, pop):
    '''
    Assigns the 
    '''
    coords = loos.Atom.coords(atom)

    dx = coords[0] + shift[0] - center[0]
    dy = coords[1] + shift[1] - center[1]
    
    # assign to bin
    if dx*dx + dy*dy <= radius**2:
        pop[0] += 1


def pbc_images(center, r, size):
    '''
    This is very brute-force at the moment. Need to figure out a more elegant 
    way of doing this.
    '''
    images = [[0,0]]
    if center[0] + r > size[0]/2:
        images.append([1,0])
    if center[1] + r > size[1]/2:
        images.append([0,1])
    if center[0] - r < -size[0]/2:
        images.append([-1,0])
    if center[1] - r < -size[1]/2:
        images.append([0,-1])
    if center[0] + r > size[0]/2 and center[1] + r > size[1]/2:
        images.append([1,1])
    if center[0] + r > size[0]/2 and center[1] + r < -size[1]/2:
        images.append([1,-1])
    if center[0] - r < -size[0]/2 and center[1] + r > size[1]/2:
        images.append([-1,1])
    if center[0] - r < -size[0]/2 and center[1] + r < -size[1]/2:
        images.append([-1,-1])
    
    return np.array(images)


def mole_fraction(pop, num_bins, species):
    '''
    Calculates mole fraction of specified species in each bin of the n*n grid 
    from a species*n*n population array. Returns an n*n array of mole fraction 
    values.
    If, for a chosen grid size, there are bins with no lipids, returns a 
    message and the mole fraction for the specific bin is set to NaN.
    '''
    molfrac = np.zeros((num_bins,num_bins))
    for i in range(num_bins):
        for j in range(num_bins):
            
            total_ij = pop[0][i][j] + pop[1][i][j] + pop[2][i][j]
            
            if total_ij == 0:
                print("Gap in membrane at bin ("+str(i)+","+str(j)+") for an "
                      +str(num_bins)+"x"+str(num_bins)+" grid.")
                molfrac[i][j] = np.nan
            else:
                molfrac[i][j] = pop[species][i][j]/total_ij
    
    return molfrac


def concentration(pop, bin_size):
    '''
    Calculates and returns a n*n concentration array from an n*n population
    array by dividing population of each bin by bin area.
    '''
    bin_area = bin_size[0]*bin_size[1] 
    
    return pop/bin_area


def enrichment(pop, bin_size, ideal_conc):
    '''
    Calculates and returns a n*n enrichment array from an n*n population
    array by dividing population of each bin by bin area and ideal 2D gas 
    concentration.
    '''
    bin_area = bin_size[0]*bin_size[1]
    
    return pop/(bin_area*ideal_conc)


def thickness(lipid_particles):
    '''
    Takes input of an array of individual lipid bead Atomic Groups and outputs
    an array of leaflet thickness at each lipid. The bead is typically the PO4
    bead below the headgroup.
    The thickness at the lipid is simply the z-distance from the membrane 
    centroid z-coordinate (which is 0 in the input trajectory).
    '''
    thickness = []
    
    for particle in lipid_particles:
        thickness.append(abs(particle.centroid().z()))
    
    return thickness


def histogram(array, hist_bins, range_min, range_max):
    '''
    Takes input of an array of real numbers of arbitrary dimensions, number of 
    histogram bins, and binning range. Calculates a numpy histogram of the
    array with specified parameters, and returns a bins*2 sized array of bin
    heights in the 0th column and bin centers in the 1st.
    '''
    hist_bin_size = (range_max - range_min)/hist_bins
    h, b = np.histogram(array, bins = np.arange(range_min, 
                                                range_max + 0.1*hist_bin_size, 
                                                hist_bin_size), density = True)
    
    # calculate bin centers
    b = b[:-1] + 0.5*(b[1]-b[0])
    
    return np.column_stack((h,b))