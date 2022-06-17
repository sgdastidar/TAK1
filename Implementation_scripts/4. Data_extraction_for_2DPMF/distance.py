import glob
import re
import os
import mdtraj as md
import numpy as np
import pandas as pd

#### Create a list of dcd files to be used for dcd calculations ####

list = []
for file in sorted(glob.glob('with*.h5')):     # Only two files are required for this calulation (withtab1.h5, withoutab1.h5)
    if re.match(r'.*ca.*', file):              # Exclude files from list having ca in the filename (applies for any pattern)
        continue                               # Breaks the loop if the 'if' condition is met
    else:
        list.append(file)                      # Else append the filename to the list 
#print(list) 


#### Read md files in small chunks for each distance calculation and append into a single array ####

for traj in list:
  DS1 = []
  DS2 = []

  #### Atom selection for distance calculation ####

  frame = md.load_frame(traj, 0)                                  # Load the first frame of the traj to get the atom indices in i, j, k etc.          
  i = frame.topology.select("protein and resSeq 154 and name CA")
  j = frame.topology.select("protein and resSeq 161 and name CA")
  k = frame.topology.select("protein and resSeq 176 and name CZ")

  d1 = frame.topology.select_pairs(k,i)                           # Define pair for computing F176(CZ)-H154(CA) distance (DS1)
  d2 = frame.topology.select_pairs(k,j)                           # Define pair for computing F176(CZ)-H154(CA) distance (DS1)

  for chunk in md.iterload(traj, chunk=1000):                     # **specify CHUNK SIZE (CHANGEABLE)** #
      #print (chunk)
      #print (chunk.time)   
   
  #### Chunk modification for analysis ####
      t = chunk[::1]                                              # **stride (CHANGEABLE)** #
      #print(t)

   #### Distance calculation ####
      D1 = md.compute_distances(t, d1)
      D2 = md.compute_distances(t, d2)
      #print(DS1, DS2)
      DS1 = np.append(DS1, D1)
      DS2 = np.append(DS2, D2)
  filename = os.path.splitext(traj)[0]                            # Extract the filename from each trajectory file into a variable
  #np.save(filename+'/'+'D1-%s'%(filename), DS1)                  # Save each trajectory's distance array (in npy format) into a separate folder (create with the same name as the 'filename').
  #np.save(filename+'/'+'D2-%s'%(filename), DS2)            

  ### Concatenate the two 1D distance arrays into a composite 2D array ###

  dismat = np.vstack((DS1, DS2)).T 
  print(dismat)
  print(filename)
  np.savetxt(filename+'/'+'dis-2D.dat', dismat, fmt='%1.4f')     # Save the 2D array into a dat file into the respective folder (with the name as 'filename').

