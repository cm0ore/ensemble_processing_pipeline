# ensemble_processing_pipeline

Pipeline for going from ensembles (.pdb) --> pymol session displaying: 1. the aligned and shared water coordinates across all of the ensembles 2. the protein residues within 2A of the shared waters --> .csv with list of protein residues within 2A of shared waters and their RMSFs

## Order of Operation

0. run ensem_rmsf.py on the .pdbs in order to generate the .pml files that contain the residue rmsfs across each ensemble. 

libtbx.python ensem_rmsf.py 7KRO_updated.pdb O "pml" > 7KRO_rmsfs.pml

1. run water_ensem_rmsf.py on the .pdbs in order to generate the .pml files that contain the average water positions across each ensemble and their rmsfs

```
libtbx.python water_ensem_rmsf.py 7KRO_updated.pdb O
```

output: 7KRO_water_rmsfs.pml and 7KRO_water_rmsfs.csv

2. run ligand_series_aligner.py with each pdb and its 7KRO_water_rmsfs.pml file in order to align all of the ensembles and their waters in the same coordinate plane. 

Note: user must change
  -path to apo .pdb and its pdb_water_rmsfs.pml and its pdb_rmsfs.pml
  -variables "directory" and "pdb_list" (pdb_list should be a .csv with the name of each .pdb file). The pdbs and their water_rmsfs.pmls and "pdb_list" should all be in directory. The protein rmsfs for each pdb should be in directory/rmsfs
  
pymol: ligand_series_aligner.py 

output: pdbname_water_aligned.pdb for each pdb. This file has only the aligned waters.

3. run pdb_water_parser.py on the aligned water .pdbs in order to find the overlapping water coordinates and generate a .pml file with the coordinates of the overlapping waters. 

libtbx.python pdb_water_parser.py apo_pdbname_water_aligned.pdb csv_with_list_of_pdbname_water_aligned.pdbs

output: aligned_waters.pml and apo_only_waters.pml

4. Load these two .pml files into your existing pymol session from step 2 

5. run X.py in your pymol session order to only display the positions of the shared waters and protein residues within 2A of the shared waters. Generates a .pdb with the close protein residues.

pymol: X.py

output: close_protein_ligandbound.pdb

6. run close_protein.py on close_protein_ligandbound.pdb in order to generate a .csv with the list of protein residues within 2A of shared waters and their RMSFs and a histogram of each amino acid. 

close_protein.py 






