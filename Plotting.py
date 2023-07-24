#

import os
import numpy as np
import matplotlib.pyplot as plt

dir = os.path.dirname(__file__)


pixel_grid = 8
hist_bins = 12
n = hist_bins - 1

labels = ["DIPC", "DPPC", "Chol"]
upper = [np.loadtxt(os.path.join(dir, "CV Data", "upper_leaflet_enrichment_"+labels[i]+".txt")) for i in range(len(labels))]
lower = [np.loadtxt(os.path.join(dir, "CV Data", "lower_leaflet_enrichment_"+labels[i]+".txt")) for i in range(len(labels))]
UU = [np.loadtxt(os.path.join(dir, "Output", "UU_Test_"+labels[i]+".txt")) for i in range(len(labels))]

L1_bins_upper = upper[0][:,1]
bin_min = L1_bins_upper[0]
bin_max = L1_bins_upper[n]

frames = int(len(upper[0][:,0])/hist_bins)

L1_hist = (upper[0][:,0] + lower[0][:,0])/2
L2_hist = (upper[1][:,0] + lower[1][:,0])/2
Ch_hist = (upper[2][:,0] + lower[2][:,0])/2

L1_h = np.empty(shape=[hist_bins, 0])
L2_h = np.empty(shape=[hist_bins, 0])
Ch_h = np.empty(shape=[hist_bins, 0])

for i in range(frames):
    t = hist_bins*i
    temp = L1_hist[t:t+hist_bins]
    L1_h = np.column_stack((L1_h,temp))
    temp = L2_hist[t:t+hist_bins]
    L2_h = np.column_stack((L2_h,temp))
    temp = Ch_hist[t:t+hist_bins]
    Ch_h = np.column_stack((Ch_h,temp))

x, y = np.meshgrid(np.linspace(0, frames, frames+1), np.linspace(bin_min, bin_max, hist_bins+1))

data = [L1_h,L2_h,Ch_h]

for element, psep, label in zip(data, UU, labels):
    
    fig = plt.figure(figsize=(30, 8), dpi=300)
    ax = fig.add_subplot(1, 1, 1)
    c = ax.pcolormesh(x, y, element, cmap='viridis')
    # set the limits of the plot to the limits of the data
    ax.axis([x.min(), x.max(), y.min(), y.max()])
    ax.set_xlabel("Frames",fontsize = 15)
    ax.set_ylabel("Enrichment",fontsize = 15)
    ax.set_title(label+" enrichment histogram over trajectory",fontsize = 20)
    fig.colorbar(c, ax=ax)
    fig.tight_layout()
    plt.savefig(os.path.join(dir,"Plots" ,"Plot_"+label+"_"+str(pixel_grid)+"_"+str(hist_bins)+".png"))
    
    
    fig = plt.figure(figsize=(30, 8), dpi=300)
    ax = fig.add_subplot(1, 1, 1)
    ax.plot(np.arange(0,frames), psep, color = 'red')
    
    # c1 = ax.vlines(psep[:], 0.0, 1.0, linestyles = 'solid', color = color, linewidth = 1.0)
    c1 = ax.axis([0, frames, -0.01, 1.01])
    ax.set_xlabel("Frames",fontsize = 15)
    # # ax.set_yticklabels(["Run3","Run2","Run1"])
    ax.set_xticklabels(np.arange(0, frames, 10))
    ax.set_yticklabels([0,1])
    ax.set_ylabel("Phase separation",fontsize = 15)
    # # ax.set_ylabel("H-bond Occupancy (binary)")
    #fig.colorbar(color)
    
    fig.tight_layout()
    plt.savefig(os.path.join(dir,"Plots" ,"UU_"+label+"_"+str(pixel_grid)+"_"+str(hist_bins)+".png"))
    
    
    # fig, ax = plt.subplots(2, 1, gridspec_kw={'height_ratios': [3, 1]}, figsize=(30, 16), dpi=300)
    
    # ax = fig.add_subplot(2, 1, 1)
    # c = ax.pcolormesh(x, y, element, cmap='viridis')
    # # set the limits of the plot to the limits of the data
    # ax.axis([x.min(), x.max(), y.min(), y.max()])
    # ax.set_xlabel("Frames",fontsize = 15)
    # ax.set_ylabel("Enrichment",fontsize = 15)
    # ax.set_title(label+" enrichment histogram over trajectory",fontsize = 20)
    # fig.colorbar(c, ax=ax)
    
    # colors = ['cyan','green']
    # color = [colors[int(sep)] for sep in psep]
    
    # ax = fig.add_subplot(2, 1, 2)
    # ax.plot(np.arange(0,frames), psep, color = 'red')
    
    # # c1 = ax.vlines(psep[:], 0.0, 1.0, linestyles = 'solid', color = color, linewidth = 1.0)
    # c1 = ax.axis([0, frames, -0.01, 1.01])
    # ax.set_xlabel("Frames",fontsize = 15)
    # # # ax.set_yticklabels(["Run3","Run2","Run1"])
    # ax.set_xticklabels(np.arange(0, frames, 10))
    # ax.set_yticklabels([0,1])
    # ax.set_ylabel("Phase separation",fontsize = 15)
    # # # ax.set_ylabel("H-bond Occupancy (binary)")
    # #fig.colorbar(color)
    
    # #fig.tight_layout()
    # # plt.show()
    # plt.savefig("Plot_"+label+"_"+str(pixel_grid)+"_"+str(hist_bins)+".png")
