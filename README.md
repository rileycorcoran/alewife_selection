# Alewife selection
Repository for the code used to process the sequencing reads from Corcoran et al. _in preparation_

All scripts in this repository start blank, and are automatically filled in with user-specific names, paths, etc. Blank scripts are located in `scripts/blank_scripts`, and the file containing all variables to be replaced is located at `scripts/wgs_pipeline_variables.sh`. This file is executable, such that, once user-specified variables are set, it will create a set of scripts corresponding to all of the user-specified information.
This is mainly useful in switching between reference genomes (either the base American shad genome, or the alternate reference based on the American shad genome created for this manuscript) without manually updating file names or paths in each script. 

The functions of all scripts can be found at `scripts/ReadMe.txt`.
