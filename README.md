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
    
