# -*- coding: utf-8 -*-
"""
Created on Fri Jul 14 00:26:49 2023

@author: saroj
"""

import os
import math
import numpy as np
from scipy import stats
from scipy.spatial import ConvexHull
import matplotlib.pyplot as plt
import argparse


dir = os.path.dirname(__file__)


def plot_ecdf(x, F, gcm, gcmf, lcm, lcmf, hullX, species, frame):
    fig = plt.figure(figsize=(20, 8), dpi=300)
    ax = fig.add_subplot(1, 1, 1)
    
    ax.plot(x, F, '.', color = 'black')
    ax.plot(gcm, gcmf, 'r--', label = 'gcm')
    ax.plot(lcm, lcmf, 'g--', label = 'lcm')
    ax.axis([min(hullX), max(hullX), 0.0, 1.0])
    ax.set_xlabel("X",fontsize = 15)
    ax.set_ylabel("CDF(X)",fontsize = 15)
    ax.set_title("GCM and LCM curves",fontsize = 20)
    
    plt.legend()
    plt.savefig(os.path.join(dir,"Plots","ecdf","Plot_gcm_lcm_"+str(frame)+"_"+species+".png"))
    plt.show()
    plt.close(fig)


def ecdf(X):
    """
    Takes a dataset of random variables, X = {x_1, ..., x_N}, sorted in 
    increasing order and returns the empirical cumulative distribution 
    function (ecdf) of X.
    ecdf(X) is given by:
        F(x_i) = (number of data points <= x_i in X)/N
        where x_i are the data points in X
    """
    N = len(X)
    F = []
    
    for point in X:
        count = N - X[::-1].index(point)
        f = count/N
        F.append(f)
    
    return F, X


def intersect(X,P0,P1):
    """
    Takes in a dataset X of random variables and two bounds, and outputs the 
    portion of the dataset within those bounds
    """
    X1 = [value for value in X if value >= P0]                  # X1 = X(x>=P0)
    X2 = [value for value in X if value <= P1]                  # X2 = X(x<=P1)
    XX = [value for value in X1 if value in X2]            # X1 intersection X2
    
    return XX


def convexhull(x, F):
    """
    Takes in ordered x and y values of an ecdf and calculates its convex hull
    points using scipy.spatial.ConvexHull. Returns the indices of the hull 
    points in the input values lists in counterclockwise order.
    """
    points = []
    for i in range(len(x)):
        points.append([x[i],F[i]])
    
    hull = ConvexHull(points)
    
    return hull.vertices


def KS(X):
    """
    Kolmogorov Smirnov uniformity test of a dataset of random variables X. 
    Returns 1 if dataset is uniform, 0 if it is non-uniform.
    """
    if len(list(np.unique(np.array(X)))) <= 2:
        return 1
    
    elif len(X)>1:
        # make a uniform distribution cdf
        low, high = min(X), max(X)
        
        # Check to avoid KS test failure for very few extremely close points.
        # If not used, it can break the scipy convex hull function due to 
        # underflow compared to its precision limit
        if max(X)-min(X) <= 0.000001:
            return 1
        
        x = np.random.uniform(low, high, 10000)
        
        # perform KS test
        res = stats.kstest(X, x)
        
        # interpret p-value output
        if res.pvalue >= 0.01:
            return 1
        else:
            return 0
    
    else:
        return 1


def gcmlcm(X, species, frame, recur):
    """
    Takes in a dataset of random variables X and calculates the ecdf of X. It
    then uses the convex hull points of the ecdf to calculate the greatest 
    convex minorant and least concave majorant of the ecdf.
    Returns lists of the gcm and lcm points.
    """
    F, x = ecdf(X)
    F.pop(0)
    x.pop(0)
    
    x_unique = np.unique(np.array(x))
    if len(x_unique) < 3:
        gcm, gcmf = [x_unique[0]], [min(F)]
        lcm, lcmf = [x_unique[1]], [max(F)]
        
        return gcm, lcm
    
    # get indices of convex hull points, create an ordered list of hull points
    indices = convexhull(x, F)
    indices = np.roll(indices, -np.argmin(indices))
    
    hullX, hullF = [], []
    for i in range(len(indices)):
        hullX.append(x[indices[i]])
        hullF.append(F[indices[i]])
    
    gcm, lcm = [], []
    gcmf, lcmf = [], []
    
    end0 = hullX.index(min(hullX))
    end1 = hullX.index(max(hullX))
    
    lcm.extend([hullX[end0], hullX[end1]])
    lcmf.extend([hullF[end0], hullF[end1]])
    for i in range(end0, end1+1):
        gcm.append(hullX[i])
        gcmf.append(hullF[i])
    for i in range(end1+1, len(indices)):
        lcm.append(hullX[i])
        lcmf.append(hullF[i])

    lcm.sort()
    lcmf.sort()
    
    # uncomment if ecdf of total data of each frame, each species is required
    # if recur == 0: plot_ecdf(x, F, gcm, gcmf, lcm, lcmf, hullX, species, frame)
    
    return gcm, lcm


def forward_search(X, PF, eL):
    
    success = False
    PF_ = [value for value in PF if value <= eL]
    eR_ind = PF.index(eL) + 2
    for i in range(eR_ind,len(PF)):
        eR = PF[i]
        XX = intersect(X, eL, eR)
        if KS(XX) == 1:
            PF_ = [eR]
            success = True
            
            return PF_, success
    
    if success == False:
        PF_ = []
    return PF_, success


def backward_search(X, PB, eR):
    
    success = False
    PB.pop()
    while len(PB) >= 1:
        XX = intersect(X, max(PB), eR)
        if KS(XX) == 1:
            PB.append(eR)
            success = True
            
            return PB, success
        
        PB.pop()
    
    return PB, success


def sufficient(P,X):
    
    success = True
    
    XX = intersect(X, P[0], P[-1])
    if KS(XX) == 1:
        P_ = [P[0],P[-1]]
        return P_, success
    
    P_ = [P[0]]
    while max(P_) != max(P):
        eR = min([value for value in P if value > max(P_)])
        
        XX = intersect(X, max(P_), eR)
        if KS(XX) == 1:
            P_.append(eR)
        
        else:
            PF, success = forward_search(X, P, max(P_))
            if success == True: P_.extend(PF)
            else:
                PB, success = backward_search(X, P_, eR)
                if success == False:
                    P_ = []
                    return P_, success
                P_ = PB
    
    return P_, success


def make_consistent(gcm, lcm, C, ind):
    """
    Takes gcm and lcm sets that lead to an inconsistent gcm U lcm set and
    makes two consistent subsets from them.
    For the first, all gcm points after the first lcm point are removed.
    For the second, all lcm points before the last gcm point are removed.
    These are returned in the list of lists C, where C[0] is the first
    and C[1] is the second consistent subset of gcm U lcm.
    The list of lists 'ind' simply contains ind[0] and ind[1] corresponding to 
    C[0] and C[1] such that gcm points in C[i] correspond to a 0 in ind[i] and
    lcm points in C[i] correspond to a 1 in ind[i].
    """
    
    temp = [value for value in gcm if value < lcm[0]]
    length = len(temp)
    C.append(temp)
    C[0].extend(lcm)
    ind.append(list(np.concatenate((np.zeros(length), np.ones(len(lcm))))))
    
    temp = [value for value in lcm if value > gcm[-1]]
    C.append(gcm)
    length = len(gcm)
    C[1].extend(temp)
    ind.append(list(np.concatenate((np.zeros(length), np.ones(len(temp))))))


def consistent(gcm, lcm):
    """
    Takes gcm and lcm lists and returns 1 or 2 possible consistent subsets of
    gcm U lcm.
    A subset of gcm U lcm is consistent when max(gcm) < min(lcm)
    """
    
    C, ind = [], []
    
    if len(gcm) == 0:
        C.append(lcm)
        ind.append(list(np.ones(len(lcm))))
    
    elif len(lcm) == 0:
        C.append(gcm)
        ind.append(list(np.zeros(len(gcm))))
    
    elif gcm[-1] < lcm[0]:
        C.append(gcm)
        length = len(gcm)
        C[0].extend(lcm)
        ind.append(list(np.concatenate((np.zeros(length), np.ones(len(lcm))))))
    
    else:
        make_consistent(gcm, lcm, C, ind)
    
    return C, ind


def UU(SG, PI, SL, X, species, frame, recur):
    
    SG_, PI_, SL_, success = SG, [], SL, False
    
    XX = intersect(X, PI[0], PI[1])
    
    # Kolmogorov-Smirnov uniformity test
    if KS(XX) == 1:
        PI_ = PI
        success = True
        return SG_, PI_, SL_, success
    
    # compute gcm and lcm subsets of X
    gcm, lcm = gcmlcm(XX, species, frame, recur)
    
    # remove last element of gcm (common with last element of lcm)
    # remove first elelent of lcm (common with first elelent of gcm)
    if len(gcm) > 1:
        gcm.pop()
    if len(lcm) > 1:
        lcm.pop(0)
    
    C, ind = consistent(gcm, lcm)
    
    for i in range(len(C)):
        pos = []
        if 1.0 in ind[i]: pos.append(ind[i].index(1.0))
        
        if len(pos) != 0:
            if len(pos) != 0 and pos[0] == 0:
                PI_ = []
                PL = C[i]
                
                PL_, success = sufficient(PL,X)
                if success == False: continue
                
                SL_.extend(PL_)
            
            elif len(pos) != 0 and pos[0] == 1:
                PI_ = C[i][0:pos[0]+1]
                PL = C[i][pos[0]+1:len(C[i])]

                if len(PL) != 0:
                    PL_, success = sufficient(PL,X)
                if success == False: continue
                
                SL_.extend(PL_)
                
            
            elif pos[0] == len(C[i])-1:
                PG = C[i][0:pos[0]-1]
                PI_ = C[i][pos[0]-1:len(C[i])]

                if len(PG) != 0:
                    PG_, success = sufficient(PG,X)
                if success == False: continue
                
                SG_.extend(PG_)
            
            else:
                PG = C[i][0:pos[0]-1]
                PI_ = C[i][pos[0]-1:pos[0]+1]
                PL = C[i][pos[0]+1:len(C[i])]
                
                PG_, success = sufficient(PG,X)
                if success == False: continue
                
                PL_, success = sufficient(PL,X)
                if success == False: continue
                
                SG_.extend(PG_)
                SL_.extend(PL_)
        
        else:
            PI_ = []
            PG = C[i]
            
            PG_, success = sufficient(PG,X)
            if success == False: continue
            
            SG_.extend(PG_)
        
        if success == True:
            if len(PI_) != 0:
                
                SG_, PI_, SL_, success = UU(SG_,PI_,SL_,X, species, frame, recur+1)
                if success == False: continue
            return SG_, PI_, SL_, success
    
    if success == False:
        return [], [], [], False #SG_, PI_, SL_, success


def UU_Test(X, species, frame):
    
    X = list(np.around(X,6))    # rounds off input to 6 decimal places. Avoids numbers like 1.5142793000000001
    X.sort()
    
    # initialize convex, intermediate, and concave segments, the success flag 
    SG, SL = [], []
    success = True
    PI = [min(X), max(X)]
    
    # call main UU test function
    SG, PI, SL, success = UU(SG, PI, SL, X, species, frame, 0)

    return SG, PI, SL, success 


# Main

if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument('--frameNum',
                        type=int,
                        default=10,
                        help='Number of frames')
    parser.add_argument('--species',
                        type=str,
                        help='List of input files', nargs='+')
    args = parser.parse_args()
    frames = args.frameNum

    species = []
    for item in range(len(args.species)):
        species.append(args.species[item])

    data = [[] for i in range(len(species))]
    
    for i in range(frames):
        upper = [np.loadtxt(os.path.join(dir, "CV Data", "Upper_leaflet_enrichment_"+species[j]+"_"+str(i)+".txt")) for j in range(len(species))]
        lower = [np.loadtxt(os.path.join(dir, "CV Data", "Lower_leaflet_enrichment_"+species[j]+"_"+str(i)+".txt")) for j in range(len(species))]
        
        CV = []
        
        #print("\nFrame: "+str(i+1)+"\n")
        for j in range(len(species)):
            # average CV for the species over leaflets
            CV.append((upper[j] + lower[j])/2)
            
            # call UU test
            SG, PI, SL, success = UU_Test(CV[j], species[j], i)
            
            # Interpret UU test results
            #print(species[j]+": "+str(success))
            
            if success == True:
                data[j].append(0)
            elif success == False:
                data[j].append(1)
    
    for i in range(len(species)):
        np.savetxt(os.path.join(dir,"Output",'UU_Test_'+species[i]+'.txt'), data[i], fmt='%.1e')
