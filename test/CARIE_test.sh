#!/bin/sh

python ../CARIE.py --GTF ./test.gtf -i ./test_R1.bam,./test_R2.bam,./control_R1.bam,./control_R2.bam --anchor 8 --length 100 --lib unstrand --read P --type All --comparison ./comparison -o ./bam --analysis U

python ../CARIE.py --GTF ./test.gtf -i ./test_R1.sam,./test_R2.sam,./control_R1.sam,./control_R2.sam --anchor 8 --length 100 --lib unstrand --read P --type All --comparison ./comparison -o ./sam --analysis U

python ../CARIE.py --GTF ./test.gtf -i ./test_R1.sam,./test_R2.sam,./control_R1.sam,./control_R2.sam --anchor 8 --length 100 --lib unstrand --read P --type Junction --comparison ./comparison -o ./Junction --analysis U


