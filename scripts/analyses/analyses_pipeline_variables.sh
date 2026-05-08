## Setting up base directory values
user_dir="/data/nsm/velottalab/riley/manuscript_pipeline"
tmp_dir="${user_dir}/data/tmp"
mkdir -p ${tmp_dir}
scripts_dir="${user_dir}/scripts/analyses"
bin_dir="${user_dir}/bin"
data_dir="/data/nsm/velottalab/riley/manuscript_pipeline/data"

reference="/data/nsm/velottalab/public_res/genomes/fish/american_shad/GCF_018492685.1/bwa_mem2/GCF_018492685.1_fAloSap1.pri_genomic_replaced.fna"

## VCF files (post hard-filtering)
vcf_in_dir="${user_dir}/data/vcf/ashad_replaced/"
vcf_in_base="all_chrom_scaffolds_ashad_replaced_allsites_bylake-vc_gatkrecHardF_snps_maxMeanDP16_softF_rmA28Q17Q1"
vcf_pairwise_in_dir="${vcf_in_dir}/pairwise"

## Ohana variables
filtercheck_dir="${scripts_dir}/data/filtercheck"
ohana_out_dir="${data_dir}/analyses/ohana/"

## Popgen Stats
vcftools_step1="_maxMissing75"
hwe_filter="1e-10"
vcf_in_suffix="${vcftools_step1}_HWE${hwe_filter}"
final_vcf="${vcf_in_base}${vcftools_step1}_HWE${hwe_filter}"
popgen_dir="${data_dir}/popgen_stats"

## Pixy-specific
genotypegvcf_dir="${data_dir}/vcf/ashad_replaced/genotypeGVCF"
gvcf_in="all_chrom_scaffolds_ashad_replaced_allsites_bylake-vc"
filtering_out="gatkrecHardF_snps_maxMeanDP16_softF_rmA28Q17Q1${vcftools_step1}_HWE${hwe_filter}"
by_chrom_dir="${vcf_in_dir}/by_chromosomes"

## BetaSan
vcf_indiv_dir="${vcf_in_dir}/indiv_lakes"


##### ----- Batch variable replacement in the corresponding scripts ----- #####
for script_file in $( ls ${scripts_dir}/blank_scripts/*.{R,sh,submit} | sed "s#${scripts_dir}/blank_scripts/##" | sort -u ); do
   ## Fill in the scripts and place them into directories that correspond to the same file structure as in the blank_scripts directory
   out_script=$( echo ${script_file} | sed -r -e "s#\.(\w+)#_upd\.\1#" -e "s#^#${scripts_dir}/#" )

   ## Base directories
   sed -e "s#USER_DIR#${user_dir}#g" \
       -e "s#TMP_DIR#${tmp_dir}#g" \
       -e "s#SCRIPTS_DIR#${scripts_dir}#g" \
       -e "s#BIN_DIR#${bin_dir}#g" \
       -e "s#DATA_DIR#${data_dir}#g" \
       ${scripts_dir}/blank_scripts/${script_file} > ${out_script}

   ## Reference
   sed -i "s#REFERENCE#${reference}#g" ${out_script}

   ## VCF files (post hard-filtering)
   sed -i -e "s#VCF_IN_DIR#${vcf_in_dir}#g" \
       -e "s#VCF_IN_BASE#${vcf_in_base}#g" \
       -e "s#VCF_PAIRWISE_IN_DIR#${vcf_pairwise_in_DIR}#g" \
       -e "s#VCF_INDIV_DIR#${vcf_indiv_dir}#g" ${out_script}

   ## Popgen Stats
   sed -i -e "s#VCFTOOLS_STEP1#${vcftools_step1}#g" \
       -e "s#HWE_FILTER#${hwe_filter}#g" \
       -e "s#VCF_IN_SUFFIX#${vcf_in_suffix}#g" \
       -e "s#FINAL_VCF#${final_vcf}#g" \
       -e "s#POPGEN_DIR#${popgen_dir}#g" ${out_script}

   ## Pixy-specific
   sed -i -e "s#GENOTYPEGVCF_DIR#${genotypegvcf_dir}#g" \
       -e "s#GVCF_IN#${gvcf_in}#g" \
       -e "s#FILTERING_OUT#${filtering_out}#g" \
       -e "s#BY_CHROM_DIR#${by_chrom_dir}#g" ${out_script}

done
