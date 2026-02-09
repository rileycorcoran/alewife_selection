## Chromosomes List
awk '{print NR  "\t" $s}' chromosomes.tsv > numbered_chromosomes.tsv
sed -i 's#^1\t#index\tchromosome\n1\t#' numbered_chromosomes.tsv

## Scaffold Group List (1 group)
awk 'NR>1 {print $1}' scaffold_groups.tsv | uniq | \
awk 'BEGIN {printf("index\tsg\n")} {printf("%d\t%s\n", ++n, $1)}' > numbered_scaff_groups.tsv

## Scaffold List (49 scaffolds)
awk '{print $2}' scaffold_groups.tsv  | sed 1d > scaffolds_tmp1.txt     # Make a new tmp file of just the scaffold names
awk '{print NR  "\t" $s}' scaffolds_tmp1.txt > numbered_scaffolds.tsv
sed -i 's#^1\t#index\tscaffold\n1\t#' numbered_scaffolds.tsv
rm scaffolds_tmp1.txt scaffolds_tmp2.txt


## Individuals in each pop separately - for loop
for lake in Amos Bride Long Pat Quon; do
   ls HAPLOCALL_DIR/${lake}*HAPLOCALL_OUT.g.vcf.gz | sed -e "s#rg.g.vcf.gz#rg#" -e "s#HAPLOCALL_DIR/##" > number_samples_tmp.txt
   awk '{print NR  "\t" $s}' number_samples_tmp1.txt > numbered_samples-${lake}.txt
   sed -i 's#^1\t#index\tsample\n1\t#' numbered_samples-${lake}.txt
   rm number_samples_tmp.txt
done

## What the output file should look like (first 3 lines):
## index   sample
## 1       Amos_10_S108_L002_paired_ashad_sorted_marked_rg
## 2       Amos_12_S118_L002_paired_ashad_sorted_marked_rg
