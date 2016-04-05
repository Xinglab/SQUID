#!/bin/sh
#$ -cwd
#$ -M sfli001@gmail.com
#$ -l h_data=4G,h_rt=336:00:00,highp
#$ -pe  dc_* 1
#$ -V

python CARIE.py --GTF ./test/Mus_musculus.Ensembl.GRCm38.78.gtf -i ./test/total_0.bam,./test/cytosol_0.bam,./test/nucleus_0.bam --anchor 6 --length 50 --lib unstrand --read S --type all --comparison ./test/comparison -o ./test

