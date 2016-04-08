## CARIE: Computational Analysis of Retained Intron Events

Requirements
------------
1. Install Python 2.7.x and corresponding versions of NumPy and
SciPy.
2. Add the Python directory to the $PATH environment variable.

Installation:
------------
The source code can be directly called from Python.

Usage:
--------------------------------

    $ python CARIE.py --GTF ./test.gtf -i ./test_R1.bam,./test_R2.bam,./control_R1.bam,./control_R2.bam --anchor 8 --length 100 --lib unstrand --read P --type All --comparison ./comparison -o ./bam --analysis U

Required Parameters:
------------
	-i/--input:
		s1.bam/s1.sam[,s2.bam/s2.sam]. Mapping results for all of samples in bam/sam format. Different samples  are sepreated by commas
	--GTF:
		The gtf file
	
     -o *            output directory
     -r *            reference genome
     -v *            VCF directory
     --gz            flag denoting VCF files are gzipped 
     --rnaedit       flag to N-mask rna editing sites
     -e              file containing RNA editing sites, can be downloaded from RADAR
                     (http://rnaedit.com/download)
    
