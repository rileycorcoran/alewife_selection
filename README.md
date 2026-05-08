# Alewife selection
Repository for the code used to process the sequencing reads from Corcoran et al. _in preparation_

All scripts in this repository start blank, and are automatically filled in with user-specific names, paths, etc. Blank scripts are located in `scripts/blank_scripts`, and the file containing all variables to be replaced is located at `scripts/wgs_pipeline_variables.sh`. This file is executable, such that, once user-specified variables are set, it will create a set of scripts corresponding to all of the user-specified information.
This is mainly useful in switching between reference genomes (either the base American shad genome, or the alternate reference based on the American shad genome created for this manuscript) without manually updating file names or paths in each script. 

The functions of all scripts can be found at `scripts/ReadMe.txt`.

## Post-processing analysis scripts
As with the processing scripts, the analysis scripts begin blank in `scripts/analyses/blank_scripts/` so that all scripts can be batch edited with specific file names and paths. The batch changing script can be located at `scripts/analyses/analyses_pipeline_variables.sh`, and once run, will create ready-to-use scripts in folders within `scripts/analyses`. The `analyses/blank_scripts` and `analyses` directories share the same file organization, being roughly separated into 5 sections:
- Population genetics statistics (FST, dXY, pi, ROH, etc.) - `popgen_stats`
- Ohana outputs (admixture plots and selection scan results) - `ohana`
- Balancing selection-specific (BetaScan analysis) - `balancing_selection`
- SLiM runs - `slim`
- American shad RNAseq analysis - `ashad_rnaseq`
General-use extra files (such as a file containing an index of chromosome names) can be found in the base `scripts/analyses` directory, and more analysis-specific files (such as a file containing indices of 1Mb intervals over which to run sitelevel pixy) can be found in their specific folders.
