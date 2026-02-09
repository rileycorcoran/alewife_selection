1_trim.submit
   - What it does:
      - FastQC raw reads, and MultiQC the files into one readable output
      - Trim the raw reads with Trimmomatic; output is sent to /manuscript_pipeline/data/trim
      - FastQC the trimmed reads, and MultiQC the files into one readable output
   - Output:
      - Trimmed reads - /manuscript_pipeline/data/trim
      - FastQC and MultiQC files - /manuscript_pipeline/data/fastqc

2_align.submit
   - What it does:
      - Aligns raw reads to the reference of choice with bwa-mem2 (either the original American shad or SNP-replaced American shad reference)
      - NOTE: This script is a template that is meant to be autofilled to create per-lake submission scripts to decrease run time
         - Per-lake scripts will be sent to /manuscript_pipeline/scripts/${alignment}_scripts/2_align
   - Output:
      - Aligned reads per sample - /manuscript_pipeline/data/bam/${alignment}
         - NOTE: Because alignment was done twice (once to each reference), aligned reads are sent to a reference-specific folder
         - NOTE: For all steps after alignment, corresponding output files are also sent to reference-specific folders

3_markDuplicates.submit
   - What it does:
      - Mark duplicate reads for each sample with Picard
      - Collect various metrics on the marked .bam files: summary information from Picard `ValidateSamFile`, overview metrics from samtools `stats`
         - Compile the samtools `stats` with MultiQC
      - Combine samtools `stats` metrics from the two reference alignments in order to compare alignment performance
   - Output:
      - Duplicate-marked .bam files - /manuscript_pipeline/data/bam/${alignment}
      - `ValidateSamFile metrics` - /manuscript_pipeline/data/bam/${alignment}/markdup_metrics
      - samtools `stats`, per individual and alignment comparison - /manuscript_pipeline/data/bam/${alignment}/samtools_stats

4_callHaplotypes.submit
   - What it does:
      - Add Read Group information to each read, necessary for `HaplotypeCaller`
         - NOTE: The actual RG names are functionally unimportant, and are just taken from parts of the original sample names returned from sequencing
      - Index the RG-added .bam files with samtools `index` (creates a .bai file)
      - Run GATK `HaplotypeCaller` per sample
      - NOTE: This script is a template meant to be autofilled to create per-sample submission scripts to decrease run time
         - Per-sample scripts will be sent to /manuscript_pipeline/scripts/${alignment}_scripts/4_callHaplotypes
         - Takes ~3 days per sample with the SBATCH parameters included in the script
   - Output
      - RG-added .bam files per sample - /manuscript_pipeline/data/bam/${alignment}
      - Haplotype-called .g.vcf.gz files per sample - /manuscript_pipeline/data/gvcf/${alignment}

5_gdbi_chrom.submit
   - What it does:
      - Run GRBImport across each chromosome (25 in total, including the mitochondrial DNA) for each population separately
      - NOTE: This script runs as an array, one array job per chromosome, and array .out and .err files are sent to /manuscript_pipeline/scripts/${alignment}_scripts/5_gdbi
   - Output
      - Genomic database files for each chromosome, per population - /manuscript_pipeline/data/gvcf/${alignment}/<lake>
5_gdbi_scaffold.submit
   - What it does:
      -	Run GRBImport across each scaffold (49 in total) for each population separately
      - NOTE: This script runs as an array, one array job per scaffold, and array .out and .err files are sent to /manuscript_pipeline/scripts/${alignment}_scripts/5_gdbi
   - Output
      -	Genomic	database files for each	scaffold, per population - /manuscript_pipeline/data/gvcf/${alignment}/<lake>

6_1_genotypeGVCF_blank.submit
   - What it does:
      - This is a template file to run GenotypeGVCFs per chromosome or scaffold interval, per each population
      - The batch file creation and submission commands are located within this template file before any of the commands
      - NOTE: All ready-to-run files created with those submission scripts will be sent to /manuscript_pipeline/scripts/${alignment}_scripts/6_genotypeGVCF
   - Output:
      - .g.vcf files per chromosome interval and scaffold, per population - /manuscript_pipeline/data/vcf/${alignment}/genotypeGVCF

6_2_genotypeGVCF_concat_index.submit
   - What it does:
      - For each population, concatonate all the chromosome and scaffold files together
      - For each population, index the concatonated files for both GATK and bcftools
      - Finally, merge all population-specific files together into a final pre-filtering .vcf
   - Output:
      - Concatonated, per-population .g.vcf files - /manuscript_pipeline/data/vcf/${alignment}/genotypeGVCF
      - GATK (.tbi) and bcftools (.csi) indices for each per-population, concatonated file - /manuscript_pipeline/data/vcf/${alignment}/genotypeGVCF
      - All-populations, merged .vcf file - /manuscript_pipeline/data/vcf/${alignment}
