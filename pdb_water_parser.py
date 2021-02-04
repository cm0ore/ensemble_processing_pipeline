from __future__ import division
from __future__ import print_function
import iotbx.pdb
import sys
import os
import subprocess
import math
import iotbx.pdb.amino_acid_codes
aa_resnames = iotbx.pdb.amino_acid_codes.one_letter_given_three_letter

#########calculate bins for the whole ensemble

ensemble_count = {}
mean_xyz = {}
xyzs = {}
added_atoms = {}
apo_only = {}

#initialize mean_xyz with base_file coordinates
base_pdb_file = sys.argv[1]
pdb_obj = iotbx.pdb.hierarchy.input(file_name=base_pdb_file)

for model in pdb_obj.hierarchy.models():
  for chain in model.chains():
    for rg in chain.residue_groups():
      if rg.atom_groups()[0].resname in aa_resnames: continue
      for ag in rg.atom_groups():
        for atom in ag.atoms():
          atom_name = (atom.xyz)
	  mean_xyz[atom_name] = atom_name
	  xyzs[atom_name] = [atom_name]
	  ensemble_count[atom_name] = 1
	  apo_only[atom_name] = atom_name

print('initialized!')

lig_pdb_file  = open(sys.argv[2]) 	#should be a csv with a list of pdb filenames that only have waters (just like the apo pdb)
for line in lig_pdb_file:
  line = line.split()
  file_name = line[0] 
  print('processing', file_name)

  pdb_obj = iotbx.pdb.hierarchy.input(file_name=file_name)


  for model in pdb_obj.hierarchy.models():
    for chain in model.chains(): 
      for rg in chain.residue_groups():
        if rg.atom_groups()[0].resname in aa_resnames: continue
        for ag in rg.atom_groups():
          for atom in ag.atoms():
  
            atom_name = (atom.xyz)
#              print(atom_name)
	    for element in mean_xyz:
              coords = element
              distance = tuple(abs(i-j) for i, j in zip(coords, atom.xyz))
              if distance[0] > 3 or distance[1] > 3 or distance[2] > 3:
                continue
              elif element != atom_name:
                mean_xyz[element] = tuple(i + j for i,j in zip(mean_xyz[element], atom.xyz))
                added_atoms[atom_name] = atom_name
                ensemble_count[element] += 1
                xyzs[element].append(atom.xyz)

            if atom_name not in added_atoms:
              mean_xyz[atom_name] = atom.xyz
              xyzs[atom_name] = [atom.xyz]
              ensemble_count[atom_name] = 1

for atom_name in mean_xyz:
  xyz_sum = mean_xyz[atom_name] #summed coordinates
  mean_xyz[atom_name] = tuple(i / ensemble_count[atom_name] for i in xyz_sum) #average coordinates for each position 

#output 

file_basename = os.path.basename(file_name).split('_')[0]
file_1 = open('aligned_waters.pml', 'w')
file_2 = open('apo_only_waters.pml', 'w')

for element in mean_xyz:

  meanx = mean_xyz[element][0]
  meany = mean_xyz[element][1]
  meanz = mean_xyz[element][2]

  if ensemble_count[element] == 1 and element in apo_only:
    print ('pseudoatom apo_only_waters, chain=ZZ, state=-1, pos=[%.4f,%.4f,%.4f]' % \
    (meanx, meany, meanz), file=file_2)

  elif ensemble_count[element] == 1:
    continue

  print ('pseudoatom shared_waters, chain=ZZ, state=-1, pos=[%.4f,%.4f,%.4f]' % \
    (meanx, meany, meanz), file=file_1)

file_1.close()
file_2.close()
lig_pdb_file.close()
