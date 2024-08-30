#!/bin/bash



###################################################### VARIANT CALLING STEPS ####################################################################


# directories
ref="ref_genome/hg38.fna"
known_sites="gatk/Homo_sapiens_assembly38.dbsnp138.vcf"
aligned_reads="aligned_reads"
reads="/run/media/punit/data4/SEQ_DATA/gatk_build/fastq_dir"
results="results_fold"
data_folder="data_folder"




# -------------------
#  QC - Run fastqc 
# -------------------

echo " QC - Run fastqc of fastqfiles"

fastqc -t 10 ${reads}/SRR062634_1.fastq -o ${reads}/
fastqc -t 10 ${reads}/SRR062634_2.fastq -o ${reads}/
  


  # --------------------------------------
# # Map to reference using BWA-MEM
# # --------------------------------------
# 
echo " Map fastq filesto reference using BWA-MEM"
# 
# # BWA index reference 
# 
# 
# # BWA alignment
 bwa mem -t 30 -R "@RG\tID:SRR062634\tPL:ILLUMINA\tSM:SRR062634" ${ref} ${reads}/SRR062634_1.fastq.gz ${reads}/SRR062634_2.fastq.gz > ${aligned_reads}/SRR062634.paired.sam
# 
# 
# 
# 
# # -----------------------------------------
# #  Mark Duplicates and Sort - GATK4
# # -----------------------------------------
# 
 echo " Mark Duplicates and Sort - GATK4"
# 
 gatk MarkDuplicatesSpark -I ${aligned_reads}/SRR062634.paired.sam -O ${aligned_reads}/SRR062634_sorted_dedup_reads.bam
# 
# 
# 
# # ----------------------------------
# #  Base quality recalibration
# # ----------------------------------
# 
# 
 echo " Base quality recalibration"
# 
# # 1. build the model
 gatk BaseRecalibrator -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --known-sites ${known_sites} -O ${data}/recal_data.table
# 
# 
# # 2. Apply the model to adjust the base quality scores
 gatk ApplyBQSR -I ${aligned_reads}/SRR062634_sorted_dedup_reads.bam -R ${ref} --bqsr-recal-file {$data_folder}/recal_data.table -O ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam 
# 
# 
# 

# 
# # ----------------------------------------------
# Call Variants - gatk haplotype caller
# # ----------------------------------------------
# 
 echo " Call Variants - gatk haplotype caller"
# 
 gatk --java-options "-Xmx16g -XX:ParallelGCThreads=10"  HaplotypeCaller -R ${ref} -I ${aligned_reads}/SRR062634_sorted_dedup_bqsr_reads.bam -O ${results}/raw_variants.vcf
# 
# 
# 
# # extract SNPs & INDELS
# 
 gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type SNP -O ${results}/raw_snps.vcf
 gatk SelectVariants -R ${ref} -V ${results}/raw_variants.vcf --select-type INDEL -O ${results}/raw_indels.vcf
# 
# 
#SNP

gatk VariantFiltration \
	-R ${ref} \
	-V ${results}/raw_snps.vcf \
	-O ${results}/filtered_snps.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 60.0" \
	-filter-name "MQ_filter" -filter "MQ < 40.0" \
	-filter-name "SOR_filter" -filter "SOR > 4.0" \
	-filter-name "MQRankSum_filter" -filter "MQRankSum < -12.5" \
	-filter-name "ReadPosRankSum_filter" -filter "ReadPosRankSum < -8.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"


#INDEL
gatk VariantFiltration \
	-R ${ref} \
	-V ${results}/raw_indels.vcf \
	-O ${results}/filtered_indels.vcf \
	-filter-name "QD_filter" -filter "QD < 2.0" \
	-filter-name "FS_filter" -filter "FS > 200.0" \
	-filter-name "SOR_filter" -filter "SOR > 10.0" \
	-genotype-filter-expression "DP < 10" \
	-genotype-filter-name "DP_filter" \
	-genotype-filter-expression "GQ < 10" \
	-genotype-filter-name "GQ_filter"



# 
# Select Variants that PASS filters
gatk SelectVariants  --exclude-filtered -V ${results}/filtered_snps.vcf  -O ${results}/analysis-ready-snps.vcf


gatk SelectVariants --exclude-filtered  -V ${results}/filtered_indels.vcf -O ${results}/analysis-ready-indels.vcf



# to exclude variants that failed genotype filters
cat analysis-ready-snps.vcf|grep -v -E "DP_filter|GQ_filter" > analysis-ready-snps-filteredGT.vcf
cat analysis-ready-indels.vcf| grep -v -E "DP_filter|GQ_filter" > analysis-ready-indels-filteredGT.vcf


# Annotate using Funcotator
gatk Funcotator \
--variant ${results}/analysis-ready-snps-filteredGT.vcf \
--reference ${ref} \
--ref-version hg38 \
--data-sources-path 42basepairs.com/download/gs/broad-public-datasets/funcotator/funcotator_dataSources.v1.7.20200521g \
--output ${results}/analysis-ready-snps-filteredGT-functotated.vcf \
--output-file-format VCF

gatk Funcotator \
--variant ${results}/analysis-ready-indels-filteredGT.vcf \
--reference ${ref} \
--ref-version hg38 \
--data-sources-path 42basepairs.com/download/gs/broad-public-datasets/funcotator/funcotator_dataSources.v1.7.20200521g \
--output ${results}/analysis-ready-indels-filteredGT-functotated.vcf \
--output-file-format VCF