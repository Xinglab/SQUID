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
	CARIE.py --main program
	rMATS.py --a program to use the CARIE output to calculate differential RI events
	Using bam files to generate all type of RI level
		$ python CARIE.py --GTF ./test.gtf -i ./test_R1.bam,./test_R2.bam,./control_R1.bam,./control_R2.bam --anchor 8 --length 100 --lib unstrand --read P --type All
	Using sam files to generate all type of RI level
		$ python CARIE.py --GTF ./test.gtf -i ./test_R1.sam,./test_R2.sam,./control_R1.sam,./control_R2.sam --anchor 8 --length 100 --lib unstrand --read P --type All
	Using CARIR.py to generate the RI count and calculate differential RI events	
		$ python CARIE.py --GTF ./test.gtf -i ./test_R1.bam,./test_R2.bam,./control_R1.bam,./control_R2.bam --anchor 8 --length 100 --lib unstrand --read P --type Junction --comparison ./comparison -o ./bam --analysis U
	rMATS.py use the output of CARIE.py to calculate differential RI events	
		$ python rMATS.py --input ./bam/result/counts_all_Junction.txt --comparison ./comparison --output ./Junction2 --p 1 --FDR 0.05 --Diff 0.05 --analysis U --C 0.0001
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
		The types of RI level calculation used. The choices are All/Junction/JunctionIntron/5Simple/3Simple/5Complex/3Complex. The details of each calcuation are explained in the README file. If type is All, all of the six types of calcuation will be carried out,but the rMATS will not performed. Please use rMATS.py to calculate p-value and FDR. The default is All
	--comparison
		A file providing the sample pairs needed to calculate the differential RI level.The format should be column 1(name of comparions), column 2 (sample 1 order in the input file,replicates seperated by commas), column 3 (sample 2 order in the input file,replicates seperated by commas). If absent, rMATS step will be skipped
	--analysis
		Type of rMATS analysis to perform. analysisType is either P or U. P is for paired analysis and U is for unpaired analysis. Default is U

Type of Calculation:
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
Examples:
------------
	python CARIE.py --GTF ./test.gtf -i ./test_R1.bam,./test_R2.bam,./control_R1.bam,./control_R2.bam --anchor 8 --length 100 --lib unstrand --read P --type All --comparison ./comparison -o ./bam --analysis U

Output list:
------------
	result:
		All of final result files are in result folder.
	counts_all_$type.txt store the inclusion and skipping counts for all of the samples
		column 1: Intron Id representing the chromosome position, start and end.
		column 2: Gene id
		column 3: Strand
		column 4: Chromosome name
		column 5: Start coordinate
		column 6: End coordinate
		column 7: Whether this intron was annotated in the gtf file as retained intron event.
		column 8: Whether this intron was overlapped with exon, or the 5' splice site was overlapped with exon or the 3' site was overlapped with exon or whether this intron is a simple intron
		column 9: Inclusion counts for all of the samples seperated by commas
		column 10: Skipping counts for all of the samples seperated by commas
		column 11: Inclusion length
		column 12: Skipping length
		column 13: Intron inclusion level for all of the samples seperated by commas

	rMATS_Result_$comparison_$type.txt store the differential RI level calculated by rMATS
		column 1: Intron Id representing the chromosome position, start and end.
		column 2: Gene id
		column 3: Strand
		column 4: Chromosome name
		column 5: Start coordinate
		column 6: End coordinate
		column 7: Whether this intron was annotated in the gtf file as retained intron event.
		column 8: Whether this intron was overlapped with exon, or the 5' splice site was overlapped with exon or the 3' site was overlapped with exon or whether this intron is a simple intron
		column 9: Inclusion counts for all replicates  of sample 1 seperated by commas
		column 10: Skipping counts for all replicates  of sample 1 seperated by commas
		column 11: Inclusion counts for all replicates  of sample 2 seperated by commas
		column 12: Skipping counts for all replicates  of sample 2 seperated by commas
		column 13: Inclusion length
		column 14: Skipping length
		column 15: p-value for differential RI level of the two samples
		column 16: FDR for differential RI level of the two samples
		column 17: RI level for sample1, replicates seperated by commas
		column 18: RI level for sample2, replicates seperated by commas
		column 19: The difference of RI level between sample1 and sample2, which is the result of average RI level of sample1 minus the average RI level of sampel2.

	test:
		A folder contains test files to run the program

	log.CARIE: Log file for running CARIE pipeline

	gtf_files:
		A folder contains different types of gtf files to run the program. Use mouse genome as examples.
	Mus_musculus.Ensembl.GRCm38.78.gtf: the ensemble gtf files. This file should be provided by user. 
	Exon_Mus_musculus.Ensembl.GRCm38.78.gtf: the gtf file contains exons only
	Intron_Mus_musculus.Ensembl.GRCm38.78.gtf: the gtf file contains intron only
	Intron_Annotated_Mus_musculus.Ensembl.GRCm38.78.gtf: the gtf file contains the attributes whether the intron was annotated as retended introns in the original gtf files
	Intron_clean_Mus_musculus.Ensembl.GRCm38.78.gtf: the gtf file contains the attributes whether the intron/5'Junction/3'Junction was overlapped with Exon and whether the intron is a simple intron. 

	counts:
		A folder contains all of the count files
	n = number of samples
	count_all.txt: a file contains the counts for all of the introns
		column 1:Intron Id representing the chromosome position, start and end.
		column 2:Gene id
		column 3:Strand
		column 4:Chromosome name
		column 5:Start coordinate
		column 6:End coordinate    column 1: Intron Id representing the chromosome position, start and end.
		column 7: Inclusion counts at 5' splice sites for sample 1
		column 8: Skipping counts at 5' splice sites for sample 1
		column 9: Inclusion counts at 3' splice sites for sample 1
		column 10: Skipping counts at 3' splice sites for sample 1
		column 11: Skipping counts of the intron for sample 1
		column 12: counts lying in the intron for sample 1
		column 13-6*(n+1): more counts for samples 2-n

	rMATS_$comparison_$type.txt
    		The input file for running rMATS.

	rMATS_$comparison_$type folder
		The folder contains the result of rMATS output.



rMATS.py
------------
Required Parameters:
-------------
	i/--input:
		The count files generated by CARIE program or provided by the user. 
	--comparion
		A file providing the sample pairs needed to calculate the differential RI level.The format should be column 1(name of comparions), column 2 (sample 1 order in the input file,replicates seperated by commas), column 3 (sample 2 order in the input file,replicates seperated by commas). If absent, rMATS step will be skipped

Optional Parameters:
------------
	-o/--output:
		The output directory. The default is current directory
	--analysis:
		Type of rMATS analysis to perform. analysisType is either P or U. P is for paired analysis and U is for unpaired analysis. Default is U
	--p
		The number of threads used to run rMATS. The default is 1
	--FDR
		The maximum FDR  used to output the significant differential RI events. The default is 0.05
	--Diff
		The minimum difference of RI level used to output the significant differential RI events. The default is 0.05

Examples:
------------
	python rMATS.py --input ./bam/result/counts_all_Junction.txt --comparison ./comparison --output ./Junction2 --p 1 --FDR 0.05 --Diff 0.05 --analysis U

Output list:
------------
	rMATS_Result_$comparison.txt store the differential RI level calculated by rMATS
   		column 1: Intron Id representing the chromosome position, start and end.
		column 2: Gene id
		column 3: Strand
		column 4: Chromosome name
		column 5: Start coordinate
		column 6: End coordinate
		column 7: Whether this intron was annotated in the gtf file as retained intron event.
		column 8: Whether this intron was overlapped with exon, or the 5' splice site was overlapped with exon or the 3' site was overlapped with exon or whether this intron is a simple intron
		column 9: Inclusion counts for all replicates  of sample 1 seperated by commas
		column 10: Skipping counts for all replicates  of sample 1 seperated by commas
		column 11: Inclusion counts for all replicates  of sample 2 seperated by commas
		column 12: Skipping counts for all replicates  of sample 2 seperated by commas
		column 13: Inclusion length
		column 14: Skipping length
		column 15: p-value for differential RI level of the two samples
		column 16: FDR for differential RI level of the two samples
		column 17: RI level for sample1, replicates seperated by commas
		column 18: RI level for sample2, replicates seperated by commas
		column 19: The difference of RI level between sample1 and sample2, which is the result of average RI level of sample1 minus the average RI level of sampel2.

	Decrease_rMATS_Result_$comparison.txt 
		The file stores the significantly decreased RI level calculated by rMATS and have the same format as rMATS_Result_$comparison.txt

	Increase_rMATS_Result_$comparison.txt 
		The file stores the significantly increased RI level calculated by rMATS and have the same format as rMATS_Result_$comparison.txt

	rMATS_$comparison.txt  
		The input file for running rMATS.

	rMATS_$comparison folder
		A folder contains the result of rMATS output.


    
