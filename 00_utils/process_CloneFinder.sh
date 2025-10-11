#! /bin/bash
source ~/.bashrc
source /home/data/cqh/miniforge3/etc/profile.d/conda.sh
conda activate clonefinder

input_file=${1}
total_dp_thd=${2}
mutant_dp_thd=${3}
maf_thd=${4}


keywd=`basename ${input_file%.*}`
wkdir=`realpath ${input_file}| xargs dirname`

# analysis env prepare
mkdir -p ${wkdir}/process/tmp; export TMPDIR="${wkdir}/process/tmp"
cp -r /home/data/cqh/software/clonefinder ${wkdir}/process/clonefinder
cloneFinder_dir="${wkdir}/process/clonefinder/"
cp ${input_file} ${cloneFinder_dir}/input_snv
cd ${cloneFinder_dir}
mkdir ${wkdir}/res_clonefinder


# run clonefinder
sed -i s/maf_thd/${maf_thd}/ options.ini
sed -i s/total_dp_thd/${total_dp_thd}/ options.ini 
sed -i s/mutant_dp_thd/${mutant_dp_thd}/ options.ini
python3 clonefinder.py snv input_snv
mv input* ${wkdir}/res_clonefinder/
cd ${wkdir}/res_clonefinder
rename s/inputsnv/${keywd}/ input*
sed s/"#"/"> "/ ${keywd}_CloneFinder.meg| awk '$1 !~ "!"' | grep -v MEGA > ${keywd}_CloneFinder.fa
rename s/_CloneFinder// *

# clean analysis env
cd ${wkdir}
rm -r ${wkdir}/process
echo ${input_file} is okk !  

