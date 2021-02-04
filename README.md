# ensemble_processing_pipeline

Pipeline for going from ensembles (.pdb) --> pymol session displaying: 1. the aligned and shared water coordinates across all of the ensembles 2. the protein residues within 2A of the shared waters --> .csv with list of protein residues within 2A of shared waters and their RMSFs

##Order of Operation

1. run water_ensem_rmsf.py on the .pdbs in order to generate the .pml files that contain the average water positions across each ensemble and their rmsfs

libtbx.python water_ensem_rmsf.py 7KRO_updated.pdb O

output: 7KRO_water_rmsfs.pml

2. 
