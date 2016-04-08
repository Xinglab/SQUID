#!/bin/sh
python ../rMATS.py --input ./bam/result/counts_all_Junction.txt --comparison ./comparison --output ./Junction2 --p 1 --FDR 0.05 --Diff 0.05 --analysis U
