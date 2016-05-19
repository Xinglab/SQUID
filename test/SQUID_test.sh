#!/bin/sh

python ../SQUID.py --GTF ./test.gtf -i ./test_R1.bam,./test_R2.bam,./control_R1.bam,./control_R2.bam --anchor 8 --length 100 --lib unstrand --read P --Cal All --RPKM gene.txt --norm 33324688,31783182,33779359,40647996 --Clean true --c1 0.05 --c2 0.05 --c3 0.05 --p 1 --Comparison ./Comparison --analysis U -o ./bam

python ../SQUID.py --GTF ./test.gtf -i ./test_R1.sam,./test_R2.sam,./control_R1.sam,./control_R2.sam --anchor 8 --length 100 --lib unstrand --read P --Cal All --RPKM gene.txt --norm 33324688,31783182,33779359,40647996 --Clean true --c1 0.05 --c2 0.05 --c3 0.05 --p 1 --Comparison ./Comparison --analysis U -o ./sam


python ../SQUID.py --GTF ./test.gtf -i ./test_R1.bam,./test_R2.bam,./control_R1.bam,./control_R2.bam --anchor 8 --length 100 --lib first --read P --Cal All --RPKM gene.txt --norm 33324688,31783182,33779359,40647996 --Clean true --c1 0.05 --c2 0.05 --c3 0.05 --p 1 --Comparison ./Comparison --analysis U -o ./bam_first

python ../SQUID.py --GTF ./test.gtf -i ./test_R1.sam,./test_R2.sam,./control_R1.sam,./control_R2.sam --anchor 8 --length 100 --lib first --read P --Cal All --RPKM gene.txt --norm 33324688,31783182,33779359,40647996 --Clean true --c1 0.05 --c2 0.05 --c3 0.05 --p 1 --Comparison ./Comparison --analysis U -o ./sam_first

