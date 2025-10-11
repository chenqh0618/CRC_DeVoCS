#! /bin/bash
source ~/.bashrc
source /home/data/cqh/miniforge3/etc/profile.d/conda.sh
conda activate snv
wkdir=/public/data/cqh_project/crc/snv_calling/ && cd $wkdir

pileup_dir=/public/data/cqh_project/crc/wgs_process/res_pileup/
res_dir=/public/data/cqh_project/crc/wgs_process/res_varscan/
ref_fa=/public/public_data/reference/genome/human_g1k_v37/human_g1k_v37.fasta
normal_sample=CRN
sample_list=`cat /public/data/cqh_project/crc/data/sample_list`
tumor_list=`cat /public/data/cqh_project/crc/data/tumor_list`


## somatic SNV calling
parallel -j 4 " varscan somatic \
  ${pileup_dir}/${normal_sample}.pileup \
  ${pileup_dir}/${sample}.pileup \
  ${res_dir}/{} \
  --output-vcf 1 \
  --min-coverage  10 \
  --min-coverage-normal 10 \
  --min-coverage-tumor 6 \
  --min-var-freq 0.2 \
  --min-freq-for-hom 0.75 \
  --normal-purity 1 \
  --tumor-purity 1 \
  --p-value 0.95 \
  --somatic-p-value 0.05 " ::: $sample_list
  


## filter
ls ${res_dir}/| grep snp | sed s/.snp.vcf// |
  parallel -j 24 "varscan somaticFilter ${res_dir}/{}.snp.vcf \
    --min_coverage 30\
    --min-reads2 5\
    --min-strands2 1\
    --min-var-freq 0.05 \
    --p-value 0.05 \
    --indel_file /{}.indel.vcf \
    --output-file hc_somatic_raw/{}.snp.vcf " 

ls hc_somatic_raw/*vcf |
  parallel -j 24 "varscan processSomatic {} \
  --min-tumor-freq 0.05\
  --max-normal-freq  0.01 \
  --p-value  0.05 " &> processSomatic.log


ln ${wkdir}/hc_somatic_raw/*snp.Somatic.hc.vcf temp/
rename s/temp_// hc_somatic_raw/*vcf
rename s/.snp.Somatic.hc// hc_somatic_raw/*vcf



## generate hc maf
VEP_DATA=/public/data/public_data/biotrainee/bigbiosoft/vep
VEP_PATH=/home/data/cqh/miniforge3/envs/snv/bin/

ls hc_somatic_raw/| grep snp.Somatic.hc.vcf | parallel -j 24 "vcf2maf.pl --input-vcf hc_somatic_raw/{} \
  --output-maf res_maf/{.}.maf --normal-id NORMAL  --tumor-id TUMOR \
  --ref-fasta ${ref_fa} \
  --vep-data  ${VEP_DATA} \
  --vep-path  ${VEP_PATH}  \
  --ncbi-build GRCh37 \
  --cache-version 112 "
  
  
  


## merge multisample vcf
high_snp=/public/data/public_data/biotrainee/annotation/chinamap/maf_1e-4_site

ls hc_somatic_raw/*snp.Somatic.hc.vcf | parallel -j 24 "bgzip {} ; tabix {}.gz"
vcf-merge hc_somatic_raw/*.vcf.gz > hc_somatic_merge/merged_hc_somatic.vcf  
vcftools --vcf merged_hc_somatic.vcf \
  --min-alleles 2 --max-alleles 2 \
  --recode --recode-INFO-all --out filtered_hc_somatic \
  --exclude-positions ${high_snp}  




## annotate vcf
table_annovar.pl  merged_hc_somatic.vcf  \
  /public/public_data/biotrainee/bigbiosoft/ANNOVAR/humandb/ \
  -buildver hg19 \
  -protocol refGene,cytoBand,exac03,avsnp147,dbnsfp30a \
  -operation g,r,f,f,f  \
  -nastring . \
  -remove -polish  -vcfinput \
  -out annovar_merged_hc_somatic









