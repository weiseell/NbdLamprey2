#Define alias for project root directory
RUN_PATH=/mnt/research/ABC/sea_lamprey_2019

mkdir $RUN_PATH/QSTAT/VCFtools
cd $RUN_PATH/genotypes
cat pops.txt | while read -r POP
do
echo '#!/bin/sh
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -t 4:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=50G
#SBATCH -J vcftools_'$POP'
#SBATCH -o '$RUN_PATH'/QSTAT/VCFtools/VCFtools.'$POP'.o

newgrp ABC
cd '$RUN_PATH'/genotypes/'$POP'
module purge
module load GCC/7.3.0-2.30  OpenMPI/3.1.1 VCFtools/0.1.15-Perl-5.28.0

#extract all genotypes (for summary stats)
#vcftools --vcf '$POP'_populations.snps.vcf --extract-FORMAT-info GT --out '$POP'_GTs

#extract depth per genotype (for summary stats)
#vcftools --vcf '$POP'_populations.snps.vcf --geno-depth --out '$POP'_depth

#mean depth per individual
#vcftools --vcf '$POP'_populations.snps.vcf --depth --out '$POP'_indiv_depth

#generate MAF summary for each SNP 
vcftools --vcf '$POP'_populations.snps.vcf --freq2 --out '$POP'_allele

#get counts so total MAF can be calculated (not just per pop)
vcftools --vcf '$POP'_populations.snps.vcf --counts2 --out '$POP'_allele
#filter so only genotypes with 8X or greater depth are included 
vcftools --vcf '$POP'_populations.snps.vcf --minDP 8 --recode --recode-INFO-all --out '$POP'_8X

#move to main file for easy access to outputs
#mv '$POP'_GTs* ../
#mv '$POP'_depth* ../
#mv '$POP'_indiv_depth* ../
mv '$POP'_8X* ../
mv '$POP'_allele* ../

scontrol show job ${SLURM_JOB_ID}' > $RUN_PATH/SHELL/VCFTools."$POP".sh

sbatch $RUN_PATH/SHELL/VCFTools."$POP".sh

done
