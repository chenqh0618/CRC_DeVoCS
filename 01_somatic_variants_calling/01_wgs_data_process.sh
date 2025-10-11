#! /bin/bash
source ~/.bashrc
source /home/data/cqh/miniforge3/etc/profile.d/conda.sh
conda activate snv

raw_fq_dir=/public1/public_data/crc_data/raw_fastq
clean_fq_dir=/public1/public_data/crc_data/clean_fastq
ref_fa=/public/public_data/biotrainee/reference/genome/human_g1k_v37/human_g1k_v37.fasta
alignment_dir=/public1/data/cqh_data/crc/wgs_process/res_align
pileup_dir=/public1/data/cqh_data/crc/wgs_process/res_pileup
sample_list_file=/public1/public_data/crc_data/sample_list


## Raw Read Quality Control
for sample in `cat ${sample_list_file}`; do
	fastp -i ${raw_fq_dir}/${sample}_R1.fastq.gz -I ${raw_fq_dir}/${sample}_R2.fastq.gz \
	 -o ${clean_fq_dir}/${sample}_R1.fastq.gz -O ${clean_fq_dir}/${sample}_R2.fastq.gz \
	 -h ${clean_fq_dir}/${sample}.fastp.html -j ${clean_fq_dir}/${sample}.fastp.json -w 64
done


## Alignment to hg19
for sample in  `cat sample_list_file` ; do
 bwa mem -t 64 -M \
    -R "@RG\tID:${sample}\tSM:${sample}\tPL:ILLUMINA\tLB:WGS\tPU:unit1" \
   ${ref_fa} \
   ${clean_fq_dir}/${sample}.R1.fastq.gz ${clean_fq_dir}/${sample}.R2.fastq.gz | samtools sort -@ 36 -o ${alignment_dir}/${sample}_sort.bam
   
  picard MarkDuplicates -I ${alignment_dir}/${sample}.bam \
	-O ${alignment_dir}/${sample}_dedup.bam \
	-M ${alignment_dir}/${sample}_dedup_metrics.txt \
	--REMOVE_DUPLICATES true
done


ls ${alignment_dir}/*_dedup.bam | parallel -j 6 "samtools index {}"


## Generate Pileup File
tumor_list_file=/public1/public_data/crc_data/tumor_list
normal_sample=CRN

for sample in  `cat tumor_list_file` ; do
	samtools mpileup -B -q 1 -f ${ref_fa} ${alignment_dir}/${normal_sample}_dedup.bam \
		${alignment_dir}/${sample}_dedup.bam > ${pileup_dir}/${sample}.pileup
		
done

