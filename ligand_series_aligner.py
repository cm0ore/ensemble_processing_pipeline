#commands for aligning structures from different ensembles

import os 
import sys

#load the apo structure and its water
cmd.load("/Users/cmoore4/Desktop/Tetrad_2020/Fraser_Rotation/Pymol/7KQO_0.6_2.5_1.pdb", "7KQO")
cmd.load("/Users/cmoore4/Desktop/Tetrad_2020/Fraser_Rotation/Pymol/7KQO_water.pml", "7KQO_water")
cmd.load("/Users/cmoore4/Desktop/Tetrad_2020/Fraser_Rotation/Pymol/7KQO_rmsfs.pml")

def aligner(pdb):
  file_basename = os.path.basename(pdb).split('_')[0] 
  pdb = directory+pdb
  cmd.load(pdb, file_basename)
  cmd.load(directory+'rmsfs/'+pdb+'_rmsf.pml')
  cmd.load("%s%s_water.pml" % (directory, file_basename), "%s_water" % file_basename)
  cmd.select("ligand", "resi 201")
  cmd.show("sticks", "ligand")
  cmd.align(file_basename, "7KQO") 
  z = cmd.get_object_matrix(file_basename)
  cmd.transform_object("%s_water" % file_basename, matrix=z)
  new_dir = directory+'water_aligned/'
  cmd.save(new_dir+"%s_water_aligned.pdb" % file_basename, "%s_water" %file_basename)


directory = '/Users/cmoore4/Desktop/Tetrad_2020/Fraser_Rotation/Pymol/best_ensembles/'  #where all of the pdbs are 
pdb_list = open(directory+'best_ensemble_filenames.txt')                                       #a list of all your ligand bound pdbs
for line in pdb_list:
  line = line.split()
  pdb = line[0]
  print(pdb)
  aligner(pdb) 

