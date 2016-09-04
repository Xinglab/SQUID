#!/bin/sh


python ../SQUID.py --GTF ./test.gtf -i ./test_R1.sam,./test_R2.sam,./control_R1.sam,./control_R2.sam --anchor 8 --length 100 --lib first --read P --Cal All --RPKM gene.txt --c1 0.05  --p 1 --Comparison ./Comparison --analysis U -o ./sam_first


python ../SQUID.py --GTF ./test.gtf -i ./test_R1.bam,./test_R2.bam,./control_R1.bam,./control_R2.bam --fasta ./test_R1_1.fq:./test_R1_2.fq,./test_R2_1.fq:./test_R2_2.fq,./control_R1_1.fq:./control_R1_2.fq,./control_R2_1.fq:./control_R2_2.fq --index=/u/home/s/shaofang/project-yxing/genome/kallisto/hg19_Ensemble74 --anchor 8 --length 100 --lib unstrand --read P --Cal All  --c1 0.05  --p 1 --Comparison ./Comparison --analysis U -o ./bam
