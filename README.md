# UU-Test
Unimodal Uniform Test implementation for LOOS


For Ashlin:

The code is run in three steps:\n
python3 membrane_enrichment.py 1DIPC.psf 1DIPC_8us.dcd 'resname=="DIPC"' 'resname=="DPPC"' 'resname=="CHOL"' 8 12

This is the membrane enrichment calculation on an 8x8 pixel grid that is used to output the raw enrichment data from each frame for each species in each leaflet.
The files will be output in the CV Data folder.

python3 UUTest.py

This reads the files from the CV Data folder and prints files for each species containing 0 or 1 for each frame in the Output folder. If you uncomment the last line in the function gcmlcm(), it will also output plots for the ecdf of each species in each frame in the folder /Plots/ecdf.

Finally, the Plotting.py file reads the files in the Output folder and plots the figures in /Plots
