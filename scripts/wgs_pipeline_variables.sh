## Setting up base directory values
user_dir="/data/nsm/velottalab/riley/manuscript_pipeline"
tmp_dir="${user_dir}/data/tmp"
mkdir -p ${tmp_dir}
scripts_dir="${user_dir}/scripts"


### 1 - Trimming & FastQC ###
raw_reads_dir="/data/nsm/velottalab/rawdata/fish/alewife_wgs_n110"
trim_dir="${user_dir}/data/trim"
mkdir -p ${trim_dir}
fastqc_dir="${user_dir}/data/fastqc"
mkdir -p ${fastqc_dir}


### 2 - .bams steps - Aligning and marking duplicates ###
alignment="ashad"
bams_dir="${user_dir}/data/bam/${alignment}"
mkdir -p ${bams_dir}

## All scripts after this are repeated for each alignment type, so scripts, .out, and .err files are put into their own directories
aligntype_scriptsDir="${scripts_dir}/${alignment}_scripts"
align_scriptsDir="${aligntype_scriptsDir}/2_align"
mkdir -p ${align_scriptsDir}

## The data was created via two separate alignments
## This loop will assign the appropriate reference file depending on the ${alignment} variable set
base_reference="/data/nsm/velottalab/public_res/genomes/fish/american_shad/GCF_018492685.1/bwa_mem2/GCF_018492685.1_fAloSap1.pri_genomic"
if [ ${alignment} == "ashad" ]; then
   ## Original ashad reference
   reference=${base_reference}
elif [ ${alignment} == "ashad_replaced" ]; then
   ## Ashad reference with fixed ashad-alewife "SNPs" replaced
   reference="${base_reference}_replaced"
fi

## Assuming the file name begins with some shared value:
align_out="_sorted"

## 3 - Mark duplicates variables
markdup_metrics_dir="${bams_dir}/markdup_metrics"
mkdir -p ${markdup_metrics_dir}
markdup_out="${align_out}_marked"
samtools_stats_dir="${bams_dir}/samtools_stats"
mkdir -p ${samtools_stats_dir}/multiqc


### 4 - Calling Haplotypes ###
haplocall_scriptsDir="${aligntype_scriptsDir}/4_haplocall"
haplocall_dir="${user_dir}/data/gvcf/${alignment}"
mkdir -p ${haplocall_dir}
haplocall_out="${markdup_out}_rg"


### 5 - GDBI variables ###
gdbi_scriptsDir="${aligntype_scriptsDir}/5_gdbi"
mkdir -p ${gdbi_scriptsDir}
gdbi_dir="${user_dir}/data/gvcf/${alignment}/gdbi"
for lake in amos bride long pat quon; do
   mkdir -p ${gdbi_dir}/${lake}
done

### 6 - GenotypeGVCF variables ###
genotypegvcf_dir="${user_dir}/data/vcf/${alignment}/genotypeGVCF"
mkdir -p ${genotypegvcf_dir}
genotypegvcf_scriptsDir="${aligntype_scriptsDir}/6_genotypeGVCF"
genotypegvcf_out="all_chrom_scaffolds_${alignment}"
genotypegvcf_merge="${genotypegvcf_out}_allsites_bylake-vc"


### 7 - General filtering variables ###
final_vcf_dir="${user_dir}/data/vcf/${alignment}"
mkdir -p ${final_vcf_dir}
out_filters="gatkrecHardF_snps_maxMeanDP16_softF_rmA28Q17Q1_maxMissing75"
filtering_out="${genotypegvcf_merge}_${out_filters}"


##### ----- Batch variable replacement in the corresponding scripts ----- #####
for script_file in $( ls ${scripts_dir}/blank_scripts/*.{R,sh,submit} | sed "s#${scripts_dir}/blank_scripts/##" | sort -u ); do
   out_script=$( echo ${script_file} | sed -r -e "s#\.(\w+)#_upd\.\1#" -e "s#^#${scripts_dir}/#" )

   ## Base directories
   sed -e "s#USER_DIR#${user_dir}#g" \
       -e "s#TMP_DIR#${tmp_dir}#g" \
       -e "s#SCRIPTS_DIR#${scripts_dir}#g" \
       ${scripts_dir}/blank_scripts/${script_file} > ${out_script}

   ## 1 - Trimming & FastQC
   sed -i -e "s#RAW_READS_DIR#${raw_reads_dir}#g" \
          -e "s#TRIM_DIR#${trim_dir}#g" \
          -e "s#FASTQC_DIR#${fastqc_dir}#g" \
          ${out_script}

   ## 2 - .bam creation & editing
   sed -i -e "s#ALIGNMENT#${alignment}#g" \
          -e "s#BAMS_DIR#${bams_dir}#g" \
          -e "s#ALIGNTYPE_SCRIPTSDIR#${aligntype_scriptsDir}#g" \
          -e "s#ALIGN_SCRIPTSDIR#${align_scriptsDir}#g" \
          -e "s#REFERENCE#${reference}#g" \
          -e "s#BASE_REF#${base_reference}#g" \
          -e "s#ALIGN_OUT#${align_out}#g" \
          ${out_script}

   ## 3 - Mark duplicates
   sed -i -e "s#MARKDUP_METRICS_DIR#${markdup_metrics_dir}#g" \
          -e "s#MARKDUP_OUT#${markdup_out}#g" \
          -e "s#SAMTOOLS_STATS_DIR#${samtools_stats_dir}#g" \
          ${out_script}

   ## 4 - Calling haplotypes
   sed -i -e "s#HAPLOCALL_SCRIPTSDIR#${haplocall_scriptsDir}#g" \
          -e "s#HAPLOCALL_DIR#${haplocall_dir}#g" \
          -e "s#HAPLOCALL_OUT#${haplocall_out}#g" \
          ${out_script}

   ## 5 - GDBI
   sed -i -e "s#GDBI_SCRIPTSDIR#${gdbi_scriptsDir}#g" \
          -e "s#GDBI_DIR#${gdbi_dir}#g" \
          ${out_script}

   ## 6 - GenotypeGVCF
   sed -i -e "s#GENOTYPEGVCF_SCRIPTSDIR#${genotypegvcf_scriptsDir}#g" \
          -e "s#GENOTYPEGVCF_DIR#${genotypegvcf_dir}#g" \
          -e "s#GENOTYPEGVCF_OUT#${genotypegvcf_out}#g" \
          -e "s#GENOTYPEGVCF_MERGE#${genotypegvcf_merge}#g" \
          ${out_script}

   ## 7- Filtering
   sed -i -e "s#FINAL_VCF_DIR#${final_vcf_dir}#g" \
          -e "s#OUT_FILTERS#${out_filters}#g" \
          -e "s#FILTERING_OUT#${filtering_out}#g" \
          ${out_script}

# NEED TO TEST AND MAKE SURE THIS PARTICULAR FILE WORKS
done

cp wgs_pipeline_variables.sh ${scripts_dir}/wgs_pipeline_variables_${alignment}.sh
