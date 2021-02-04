#takes a pdb and its rmsfs 
#I'm inputing the close_protein_apo_only.pdb or the close_protein_ligandbound.pdb
#optionally takes the original pdb as a third argument to tell you which residues aren't close to water


from __future__ import print_function
import iotbx.pdb
from iotbx.pdb import hierarchy
import sys
import os
import math
import iotbx.pdb.amino_acid_codes

aa_resnames = iotbx.pdb.amino_acid_codes.one_letter_given_three_letter

file_name = sys.argv[1]
protein_rmsf_file = open(sys.argv[2])
rmsfs = {}
for line in protein_rmsf_file:
  line = line.split(",")
  if line[0] == 'dataset':
    continue
  resnum = int(line[2])
  rmsf = line[3]
  rmsfs[resnum] = rmsf

protein_rmsf_file.close()
#print(rmsfs)

pdb_obj = iotbx.pdb.hierarchy.input(file_name=file_name)
close_resnums = []
close_protein = {}

for model in pdb_obj.hierarchy.models():
  for chain in model.chains():
    for rg in chain.residue_groups():
      if rg.atom_groups()[0].resname in aa_resnames: 
        for ag in rg.atom_groups():
          resid = (ag.resname, int(rg.resseq))

          if not close_resnums:
            close_resnums = [resid[1]]
             
          if resid not in close_protein:         
            close_protein[resid] = rmsfs[resid[1]]
            close_resnums.append(resid[1])

#file_basename = os.path.basename(file_name).split('_')[0] 
outfile = open('protein_close_to_water.csv', 'w')     

print('residue within 2A of aligned waters, rmsf', file=outfile) 
for residue in sorted(close_protein, key=lambda x: x[1]):
  print(residue, close_protein[residue], file=outfile)

if len(sys.argv) > 3:
  original_pdb = (sys.argv[3])
  pdb_obj = iotbx.pdb.hierarchy.input(file_name=original_pdb)
  
  print('residues NOT within 2A of aligned waters', file=outfile)
  all_protein = {}
  for model in pdb_obj.hierarchy.models():
    if int(model.id) == 1:
      for chain in model.chains(): 
        if chain.id == "A":
          for rg in chain.residue_groups():
            for ag in rg.atom_groups():
              resid = (ag.resname, int(rg.resseq))    
	      all_protein[resid] = rmsfs[resid[1]]

  for residue in sorted(all_protein, key=lambda x: x[1]):
    if residue not in close_protein:
      print(residue, all_protein[residue], file=outfile)

outfile.close()
