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
Optional Parameters:
------------	
	--o/--output:
		The output directory. The default is current directory
	--lib
		The library type with choices of unstrand/first/second. The details are explained at the parameter of library-type in tophat2. The default is unstrand
	--read
		The sequencing strategy of producing reads with choices of P/S. The default is P
	--length
		The reads length in nucleotide. The default length is 100
	--anchor
		The anchor length in nucleotide. The program will only count reads spanning junctions with at least this anchor length on each side. The default is 8
	--type
		The types of RI level calculation used. The choices are All/Junction/JunctionIntron/5Simple/3Simple/5Complex/3Complex. The details of each calcuation are explained in the README file. If type is All, all of the six types of calcuation will be carried out,but the rMATS will not performed. The default is All
	--comparison
		A file providing the sample pairs needed to calculate the differential RI level.The format should be column 1(name of comparions), column 2 (sample 1 order in the input file,replicates seperated by commas), column 3 (sample 2 order in the input file,replicates seperated by commas). If absent, rMATS step will be skipped
	--analysis
		Type of rMATS analysis to perform. analysisType is either P or U. P is for paired analysis and U is for unpaired analysis. Default is U

Type of calculation:
------------	
	Junction: 
		The skipping counts are the reads that connect the intron start and end sites. 
		The inclusion counts are the reads that either span the intron start or intron end sites.
	JunctionIntron:
		The skipping counts are the reads that connect the intron start and end sites.
                The inclusion counts are the reads that either span the intron start or intron end sites and the reads that locate with the intron region.
	5Simple: 
                The skipping counts are the reads that connect the intron start and end sites.
                The inclusion counts are the reads that span the intron start sites.
        3Simple:
                The skipping counts are the reads that connect the intron start and end sites.
                The inclusion counts are the reads that span the intron end sites.
	5Complex: 
                The skipping counts are the reads spliced at the intron start sites.
                The inclusion counts are the reads that span the intron start sites.
        3Complex:
                The skipping counts are the reads spliced at the intron  end sites.
                The inclusion counts are the reads that span the intron end sites.

    
