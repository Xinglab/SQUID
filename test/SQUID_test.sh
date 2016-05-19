#!/bin/sh

python ../SQUID.py --GTF ./test.gtf -i ./test_R1.bam,./test_R2.bam,./control_R1.bam,./control_R2.bam --anchor 8 --length 100 --lib unstrand --read P --Cal All --Comparison ./Comparison -o ./bam --analysis U --RPKM gene.txt --norm 33324688,31783182,33779359,40647996 --c1 0.0001 --c2 0.1 --p 1 --Clean true --c3 0.0001

python ../SQUID.py --GTF ./test.gtf -i ./test_R1.sam,./test_R2.sam,./control_R1.sam,./control_R2.sam --anchor 8 --length 100 --lib unstrand --read P --Cal All --Comparison ./Comparison -o ./sam --analysis U --RPKM gene.txt --norm 33324688,31783182,33779359,40647996 --c1 0.0001 --c2 0.1 --p 1 --Clean true --c3 0.0001

python ../SQUID.py --GTF ./test.gtf -i ./test_R1.bam,./test_R2.bam,./control_R1.bam,./control_R2.bam --anchor 8 --length 100 --lib first --read P --Cal All --Comparison ./Comparison -o ./bam_first --analysis U --RPKM gene.txt --norm 33324688,31783182,33779359,40647996 --c1 0.0001 --c2 0.1 --p 1 --Clean true --c3 0.0001

python ../SQUID.py --GTF ./test.gtf -i ./test_R1.sam,./test_R2.sam,./control_R1.sam,./control_R2.sam --anchor 8 --length 100 --lib first --read P --Cal All --Comparison ./Comparison -o ./sam_first --analysis U --RPKM gene.txt --norm 33324688,31783182,33779359,40647996 --c1 0.0001 --c2 0.1 --p 1 --Clean true --c3 0.0001


