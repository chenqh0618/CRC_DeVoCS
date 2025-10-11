#! /bin/bash
source ~/.bashrc
source /home/data/cqh/miniforge3/etc/profile.d/conda.sh
conda activate cnv

res_dir=/public/data/cqh_project/crc/wgs_process/res_seqz/
pileup_dir=/public/data/cqh_project/crc/wgs_process/res_pileup/
ref_fa=/public/public_data/reference/genome/human_g1k_v37/human_g1k_v37.fasta
normal_sample=CRN
sample_list_file=/public/data/cqh_project/crc/data/sample_list



# process a fasta file to produce a gc wiggle track file
gc_wiggle="/public/data/public_data/biotrainee/reference/genome/human_g1k_v37/human_g1k_v37_gc50based.wig.gz"
sequenza-utils gc_wiggle \
  -f ${ref_fasta} \
  -w 50 \
  -o ${gc_wiggle}


# process pileup and wiggle files to produce a seqz file 
parallel -j 5 "sequenza-utils bam2seqz --pileup \
  --normal ${pileup_dir}/${normal_sample}.pileup.gz \
  --tumor  ${pileup_dir}/{1}.pileup.gz \
  -gc ${gc_wiggle} \
  --output ${res_dir}/raw_{1}.seqz.gz \
  --hom 0.9 \
  --het 0.25 \
  --qlimit 20 \
  --qformat sanger \
  -N 20 " ::: `cat ${sample_list_file}`  

parallel -j 5 "sequenza-utils seqz_binning \
  -s ${res_dir}/{} \
  --window 1000000 \
  -o ${res_dir}/binned_{} " ::: `ls ${res_dir} | grep seqz.gz | grep -v raw`





