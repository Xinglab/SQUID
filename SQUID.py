#!/bin/python
import getopt,copy,re,os,sys,logging,time,datetime;
options, args = getopt.getopt(sys.argv[1:], 'i:o:',['input=','GTF=','fasta=','index=','output=','lib=','read=','length=','anchor=','Cal=','RPKM=','Comparison=','analysis=','c1=','p='])
input ='';
GTF ='';
fasta='';
index=''
output ='.'
lib ='unstrand'
read ='P'
length = 100
anchor = 8
Cal ='All'
RPKM=''
Comparison = ''
analysis = 'U'
c1 = 0.0001
p = 1
for opt, arg in options:
	if opt in ('-i','--input'):
		input = arg
	elif opt in ('--GTF'):
               	GTF = arg
	elif opt in ('--fasta'):
                fasta = arg
	elif opt in ('--index'):
                index = arg
	elif opt in ('-o','--output'):
		output = arg
	elif opt in ('--lib'):
		lib = arg
	elif opt in ('--read'):
                read = arg
	elif opt in ('--length'):
                length = int(arg)
	elif opt in ('--anchor'):
                anchor = int(arg)
	elif opt in ('--Cal'):
                Cal = arg
	elif opt in ('--RPKM'):
                RPKM = arg
	elif opt in ('--Comparison'):
               	Comparison = arg
        elif opt in ('--analysis'):
                analysis = arg
	elif opt in ('--c1'):
                c1 = float(arg)
	elif opt in ('--p'):
                p = int(arg)
if (not input or not GTF):
	print "Not enough parameters!"
	print "Program : ", sys.argv[0]
	print "          A python program to calculate the retained intron level and differential retained introns.\n"
	print "Usage :", sys.argv[0], " -i/--input: s1.bam/s1.sam[,s2.bam/s2.sam,...]. Mapping results for all of samples in bam/sam format. Different samples are sepreated by commas;"
	print "Usage :", sys.argv[0], " --GTF: The gtf file;"
	print "Usage :", sys.argv[0], " fasta: s1_1.fq[:s1_2.fq][,s1_1.fq[:s2_2.fq],...]. The raw sequencing reads in fasta or fastq format;"
	print "Usage :", sys.argv[0], " index: the path to the kallisto index that is required to run kallisto from raw reads;"
	print "Usage :", sys.argv[0], " -o/--output: The output directory. The default is current directory;"
	print "Usage :", sys.argv[0], " --lib: The library type with choices of unstrand/first/second. The details are explained in the parameter of library-type in tophat2. The default is unstrand;"
	print "Usage :", sys.argv[0], " --read: The sequencing strategy of producing reads with choices of P/S. The default is P;"
	print "Usage :", sys.argv[0], " --length: The read length in nucleotide. The default length is 100;"
	print "Usage :", sys.argv[0], " --anchor: The anchor length in nucleotide. The program will only count reads spanning junctions with at least this anchor length on each side. The default is 8;"
	print "Usage :", sys.argv[0], " --Cal: Which  part of the program user choose to run, the choices are All/count/rMATS. All means run the whole program, count means only run the PI value calculation part, rMATS means only run the differential analysis of retained intron.  The default is All;"
	print "Usage :", sys.argv[0], " --RPKM: a file providing the RPKM value for each sample, the first column is gene ID with the following column being the RPKM value for each sample. If RPKM value is empty, the run of cufflinks will be called to generate RPKM value;"
	print "Usage :", sys.argv[0], " --Comparison: A file providing the sample pairs to calculate the differential RI level.The format should be column 1(name of comparions), column 2 (sample 1 order in the input file,replicates seperated by commas), column 3 (sample 2 order in the input file,replicates seperated by commas),column 4(type of PI value use to perform rMATS), column 5 (optional, if present as 'pool', the replicates are combined together). If absent, rMATS step will be skipped;"
	print "uasge: ", sys.argv[0], " --analysis: Type of rMATS analysis to perform. analysisType is either P or U. P is for paired analysis and U is for unpaired analysis. Default is U;"
	print "Usage :", sys.argv[0], "--c1: The cutoff of splicing difference using Junction method. The cutoff used in the null hypothesis test for differential splicing. The default is 0.0001;"
        print "Usage :", sys.argv[0], " --p: The number of threads used to run rMATS. The default is 1;"
	print datetime.datetime.now()
	print "Author  : Shaofang Li"
	print "Contact : sfli001@gmail.com"
	sys.exit()

if (not os.path.exists(output)):
        os.system("mkdir %s" % output)

### setting up the logging format 
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=output+'/log.SQUID' + str(datetime.datetime.now())+'.txt' ,
                    filemode='w')

def listToString(ss):
  Str = '';
  for a in ss:
    Str += a+' ';
  return Str;

##### Getting Start Time ######
logging.debug('Start the program with [%s]\n', listToString(sys.argv));
startTime = time.time();

#get path of the main program
path = os.path.abspath(os.path.dirname(__file__));
##get the path  of the python  programs
bin_path = "%s/bin" % path
samples = input.split(",")
num = len(samples)
##get the path  of count files
count_path = "%s/counts" % output
if (not os.path.exists(count_path)):
	os.system("mkdir %s" % count_path)
##get the path of the results
output_path = "%s/result" % output
if (not os.path.exists(output_path)):
	os.system("mkdir %s" % output_path)
normF=[0] * len(samples)
if(Cal=="All" or Cal=="count"):
	if(Cal=="All"):
		logging.debug("Run the whole program\n")
	else:
		logging.debug("Run the PI value calculation parts only\n")
	##make directory for gtf files
	gtf_path = "%s/gtf_files" % output
	if (not os.path.exists(gtf_path)):
		os.system("mkdir %s" % gtf_path)
	gtf = re.sub(".*/","",GTF)
	cmd = "cp %s %s/%s" %(GTF,gtf_path,gtf)
	os.system(cmd)
	
	##use the awk command to generate Exon.gtf file
	logging.debug("################ Generating  the intron gtf files #######################\n");
	cmd = "less %s/%s | awk '{if($3==\"exon\"){print $_}}'> %s/Exon_%s" %(gtf_path, gtf, gtf_path, gtf)
	os.system(cmd)
	logging.debug("gtf_files\Exon_" + gtf);

	##use the Intron_gtf.py to generate Intron.gtf file
	cmd = "python %s/Intron_gtf.py --gtf %s --path %s" %(bin_path, gtf,gtf_path)
	os.system(cmd)
	logging.debug("gtf_files\Intron_" + gtf);

	##Get the annotated and clean intron
	cmd = "python %s/Annotated_Intron.py --gtf %s --path %s" %(bin_path, gtf,gtf_path)
	os.system(cmd)
	logging.debug("gtf_files\Intron_Annotated_" + gtf);
	cmd = "python %s/Attri_Intron.py --gtf %s --path %s --strand %s" %(bin_path, gtf,gtf_path,lib)
	os.system(cmd)
	logging.debug("gtf_files\Intron_clean_" + gtf);
	logging.debug("#########################################################################\n");

	##generate the counts files
	if (re.search("\.bam$",input)):
		logging.debug("Using bam files to generate count files")
		cmd = "python %s/Intron_countall_strandBam.py --gtf %s/Intron_%s --length %s --anchor %s --bam %s -o %s/count_all.txt --lib %s --read %s --Total %s/Total.txt" %(bin_path,gtf_path, gtf, length, anchor, ",".join(samples), count_path, lib, read, count_path)
		logging.debug("Generate the count file")
		logging.debug(cmd)
		os.system(cmd)
	else:
		logging.debug("Using sam files to generate count files")
		cmd = "python %s/Intron_countall_strandSam.py --gtf %s/Intron_%s --length %s --anchor %s --sam %s -o %s/count_all.txt --lib %s --read %s --Total %s/Total.txt" %(bin_path,gtf_path, gtf, length, anchor, ",".join(samples), count_path, lib, read, count_path)
		logging.debug("Generate the count file")
		logging.debug(cmd)
		os.system(cmd)
	logging.debug("Finish generating the count files\n")
	logging.debug("#########################################################################\n");
	
	fr = open("%s/Total.txt" % count_path)
	info = fr.readline()
	fr.close()
	normF = map(int,info.strip().split("\t"))
	###generate the attributes of the intron
	intron_anno = dict()
	fr = open("%s/Intron_Annotated_%s" % (gtf_path, gtf))
	for info in fr:
		a = info.strip().split("\t")
		ann= re.sub('.*annotated_IR "|\".*','',a[8])
		key = "%s_%s_%s" % (a[0],a[3],a[4])
		if(key in intron_anno):
			if (ann == "true"):
				intron_anno[key] = ann
		else:
			intron_anno[key] = ann
	fr.close()
	intron_clean = dict()
	fr = open("%s/Intron_attri_%s" % (gtf_path, gtf))
	for info in fr:
		a = info.strip().split("\t")
		C1= re.sub('.*clean "|\".*','',a[8])
		cc = re.sub('.*clean_simple "|\".*','',a[8])
		key = "%s_%s_%s" % (a[0],a[3],a[4])
		if(key in intron_clean):
			if(C1=="false"):
				intron_clean[key][0]= C1
		else:
			intron_clean[key] = [C1,cc]
	fr.close()
	logging.debug("Finish generating the attributes of each introns\n")
	###generate the output 

	logging.debug("Apply Junction method")
	fr1 =open("%s/counts/count_all.txt" % output)
	##generate the counts info for all of the intron
	fw = open("%s/counts_all_Junction.txt" %output_path, "w")
	fw.write("Intron_id\tgene_id\tstrand\tchr\tstart\tend\tannotated\tclean\tinclusion_counts\tskip_counts\tIncFormLen\tSkipFormLen\tInclusionlevel\n")
	l = len(samples)	
	for info1 in fr1:	
		a1 = info1.strip().split("\t")
		if(a1[3]=="false"):
			intron_clean[a1[0]][0]= "false"
		skp = [0] * l
		inc = [0] * l
		in_level=[0] *l
		sk_l = length - 2*anchor +1
		in_l = 2* (length - 2*anchor +1 )
		for i in range(0,l):
			inc[i] = str(int(a1[i*6+7]) + int(a1[i*6+9]))
			skp[i] = a1[i*6+11]
			if(a1[i*6 +8]!= a1[i*6+11] or a1[i*6 +10]!= a1[i*6+11]):
                                intron_clean[a1[0]][1]="false"
			if((inc[i] + skp[i]) > 0) & ((float(inc[i])/in_l + float(skp[i])/sk_l) > 0):	
				in_level[i]= str((float(inc[i])/in_l)/(float(inc[i])/in_l + float(skp[i])/sk_l))
			else:
				in_level[i] ='NA'
		fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (a1[0],a1[1],a1[2],a1[4], a1[5],a1[6],intron_anno[a1[0]],",".join(intron_clean[a1[0]]), ",".join(inc), ",".join(skp), in_l,sk_l,",".join(in_level)))
	fw.close()
	fr1.close()
	logging.debug("Done generating PI_Junction\n")

	if(not os.path.exists(RPKM) or RPKM ==''):
		RPKM_path = "%s/RPKM" % output
		if (not os.path.exists(RPKM_path)):
        		os.system("mkdir %s" % RPKM_path)
		if(not os.path.exists(index) or fasta ==''):
			l_type = "fr-unstranded"
			if(lib =="first"):
				l_type = "fr-firststrand"
			if(lib =="second"):
				l_type = "fr-secondstrand"
			gene_RPKM =dict()
			
			for ss in range(0, len(samples)):
				cmd = "cufflinks --GTF %s/%s -p 1 --library-type %s --multi-read-correct -o %s/cufflinks_%s %s" %(gtf_path, gtf, l_type, RPKM_path, ss, samples[ss])
				logging.debug(cmd)
				os.system(cmd)
				fr = open("%s/cufflinks_%s/genes.fpkm_tracking" % (RPKM_path, ss))
				info = fr.readline()
				for info in fr:
					a = info.strip().split("\t")
					if(gene_RPKM.has_key(a[3])):
						gene_RPKM[a[3]][ss]=a[9]
					else:
						gene_RPKM[a[3]]=[0] * len(samples)
						gene_RPKM[a[3]][0] = a[9]
				fr.close()

			RPKM = "%s/gene_exp.txt" % RPKM_path
			fw = open(RPKM, "w")
			for rp in gene_RPKM:
				fw.write("%s\t%s\n" % (rp, "\t".join(gene_RPKM[rp])))
			fw.close()
		else:
			trans = dict()
			fr = open("%s/%s" % (gtf_path, gtf))
			for info in fr:
				a = info.strip().split("\t")
				transcript_id = re.sub('.*transcript_id "|\".*','',a[8])
				gene_id = re.sub('.*gene_id "|\".*','',a[8])
				trans[transcript_id]= gene_id
			fr.close()
			gene_RPKM =dict()
			fq = fasta.split(",")
			for ss in range(0, len(samples)):
                                cmd = "kallisto quant --index=%s --output-dir=%s/kallisto_%s --threads=4 --plaintext %s" % (index, RPKM_path, ss, re.sub(":"," ",fq[ss]))
                                logging.debug(cmd)
                                os.system(cmd)	
				fr = open("%s/kallisto_%s/abundance.tsv" % (RPKM_path,ss))
				info = fr.readline()
				for info in fr:
					a = info.strip().split("\t")
        				if(a[0] in trans):
                		#		print a[0],trans[a[0]]
						if( trans[a[0]] in gene_RPKM):
							gene_RPKM[trans[a[0]]][ss] += float(a[3])/float(a[2]) * 1000 * 1000000/ normF[ss]
						else:
							gene_RPKM[trans[a[0]]] = [0] * len(samples)
							gene_RPKM[trans[a[0]]][ss] = float(a[3])/float(a[2]) * 1000 * 1000000/ normF[ss]
				fr.close()
			RPKM = "%s/gene_exp.txt" % RPKM_path
			fw = open(RPKM, "w")
			for rp in gene_RPKM:
				fw.write("%s\t%s\n" % (rp, "\t".join(str(x) for x in gene_RPKM[rp])))
			fw.close()		
	if(os.path.exists(RPKM)):
		Gene = dict()
		fr = open(RPKM)
		info = fr.readline()
		for info in fr:
			a = info.strip().split("\t")
			Gene[a[0]] = a[1:]
		fr.close()
		intron_obs= dict()
		intron_exp = dict()
		fr = open("%s/counts/count_all.txt" % output)
		for info in fr:
			a = info.strip().split("\t")
			ri_obs = []
			ri_exp = []
			length = int(a[6])-int(a[5]) +1
#			gene_id = re.sub(",.*","",a[1])
			gene_id = a[1].split(",")
			for i in range(0,num):
				ri_obs.append(int(a[i*6 +12]) + int(a[i*6 +9]) + int(a[i*6 +7]))
				gene_exist = "false"
                                gene_expression = 0
                                
                                for ID in gene_id:
                                        if(ID in Gene):
                                                gene_exist = "true"
                                                gene_expression += float(Gene[ID][i] )
                                
                                if (gene_exist =="true"):
                                        ri_exp.append(int(gene_expression*length * normF[i] /(1000 * 1000000)))
                                else:
                                        ri_exp.append("NA")
			intron_obs[a[0]] = ri_obs
			intron_exp[a[0]] = ri_exp
		fr.close()
		fw1 =open("%s/counts/count_all_Density.txt" % output,"w")
		fr = open("%s/counts/count_all.txt" % output)
		for info in fr:
			a = info.strip().split("\t")
			length = int(a[6])-int(a[5]) +1
			fw1.write("%s\t%s\t%s\t%s\n" % (a[0],length,"\t".join(str(x) for x in intron_obs[a[0]]), "\t".join(str(x) for x in intron_exp[a[0]])))
		fw1.close()
		fr.close()
		logging.debug("Generate the Density counts")
		
		logging.debug("Apply Density method")
		fr1 =open("%s/counts/count_all_Density.txt" % output)
		fr2 =open("%s/counts_all_Junction.txt" % output_path)
		##generate the counts info for all of the intron
		fw = open("%s/counts_all_Density.txt" %output_path, "w")
		fw.write("Intron_id\tgene_id\tstrand\tchr\tstart\tend\tannotated\tclean\tobservation_counts\texp-observation_counts\tIncFormLen\tSkipFormLen\tInclusionlevel\n")
		l = len(samples)
		info1 = fr1.readline()
		info2 = fr2.readline()
		info2 = fr2.readline()
		while(info1):        
			a1 = info1.strip().split("\t")
			a2 = info2.strip().split("\t")
			if(re.search("NA",info1)):
                        	info1 = fr1.readline()
				info2 = fr2.readline()
				continue
			obs = [0] * l
			Exp_obs = [0] *l
			pi  = [0] *l
			for i in range(0,l):
				obs[i] =int( a1[i +2])
				if(int(a1[i+l+2]) < int( a1[i +2])):
					Exp_obs[i]= 0
					pi[i] =1
				else:
					Exp_obs[i] = int(a1[i+l+2]) - int( a1[i +2])
					if(int(a1[i+l+2]) >0):
						pi [i] = float(a1[i +2])/int(a1[i+l+2])
					else:
						pi[i] ='NA'
			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (a1[0],a2[1],a2[2],a2[3], a2[4],a2[5],intron_anno[a1[0]],",".join(intron_clean[a1[0]]), ",".join(str(x) for x in obs), ",".join(str(x) for x in Exp_obs), a1[1],a1[1],",".join(str(x) for x in pi)))
			info1 = fr1.readline()
                        info2 = fr2.readline()
		fw.close()
		fr1.close()
		fr2.close()
		logging.debug("Done generating PI_Density\n")

logging.debug("#########################################################################\n");

if(Cal=="count"):
	currentTime = time.time();
        runningTime = currentTime-startTime; ## in seconds
        logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60));
        logging.debug("Program finished")
        sys.exit()
if (not os.path.exists(Comparison) or not Comparison):
	logging.debug("The comparion file was not provided")
	currentTime = time.time();
        runningTime = currentTime-startTime; ## in seconds
        logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60));
	logging.debug("Program finished")
	sys.exit()
logging.debug("Run rMATS ")
if (not os.path.exists("%s/rMATS_files" % output)):
        os.system("mkdir %s/rMATS_files" % output)
fr =open(Comparison)
for info in fr:
	a = info.strip().split()
	ss1= a[1].split(",")
        ss2= a[2].split(",")
	if(a[3]=="Junction" or a[3]=="All"):
		fw = open("%s/rMATS_files/rMATS_%s_Junction.txt" % (output,a[0]),"w")	
		fw.write("ID\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\n")
		fr1 = open("%s/counts_all_Junction.txt" % (output_path))
		info1 = fr1.readline()
		for info1 in fr1:
			a1 = info1.split("\t")
			SUM = 0	
			skp1 = [0] * len(ss1)
			inc1 = [0] * len(ss1)
			skp2 = [0] * len(ss2)
			inc2 = [0] * len(ss2)
			in_level = a1[8].split(",")
			sk_level = a1[9].split(",")
			for i in range(0,len(ss1)):
				inc1[i] = in_level[int(ss1[i])-1]
				skp1[i] = sk_level[int(ss1[i])-1]
				SUM+= int(inc1[i]) + int(skp1[i])
			for i in range(0,len(ss2)):
				inc2[i] = in_level[int(ss2[i])-1]
				skp2[i] = sk_level[int(ss2[i])-1]	
				SUM+= int(inc2[i])+ int(skp2[i])
			if(SUM > 0):
				if(len(a) > 3 and a[4] =="pool"):
					fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (a1[0],sum(map(lambda x:int(x),inc1)),sum(map(lambda x:int(x),skp1)),sum(map(lambda x:int(x),inc2)),sum(map(lambda x:int(x),skp2)),a1[10],a1[11]))
				else:
					fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (a1[0],",".join(inc1),",".join(skp1), ",".join(inc2),",".join(skp2),a1[10],a1[11]))
		fr1.close()
		fw.close()
		cmd ="%s/MATS/rMATS.sh -d %s/rMATS_files/rMATS_%s_Junction.txt -o %s/rMATS_files/rMATS_%s_Junction -p %s  -t %s -c %s" %(bin_path,output,a[0],output,a[0],p,analysis, c1)
		print cmd
		os.system(cmd)
		logging.debug(cmd)
		logging.debug("Done running the rMATS for " + a[0]+ " using Junction methods")
		fr1 = open("%s/counts_all_Junction.txt" % (output_path))
		fr2 = open("%s/rMATS_files/rMATS_%s_Junction/rMATS_Result.txt" %(output,a[0]))
		fw = open("%s/result/rMATS_Result_%s_Junction.txt" %(output, a[0]),"w")
		fw.write("Intron_id\tgene_id\tstrand\tchr\tstart\tend\tannotated\tclean\tinclusions_counts_SAMPLE1\tskip_counts_SAMPLE1\tinclusions_counts_SAMPLE2\tskip_counts_SAMPLE2\tInclusion_length\tSkipping_length\tPValue\tFDR\tIncLevel_SAMPLE1\tIncLevel_SAMPLE2\tIncLevelDifference\n")
                info1 = fr1.readline()
                info2 = fr2.readline()
                RE = dict()
                for info2 in fr2:
                        a2 = info2.strip().split("\t")
                        RE[a2[0]] = "\t".join(a2[1:])
                for info1 in fr1:
                        a1 = info1.strip().split("\t")
                        if(RE.has_key(a1[0])):
                                RE[a1[0]] = "\t".join(a1[0:8]) + "\t" + RE[a1[0]]
                fr1.close()
                fr2.close()
                fr2 = open("%s/rMATS_files/rMATS_%s_Junction/rMATS_Result.txt" %(output,a[0]))
                info2 = fr2.readline()
                for info2 in fr2:
                        a2 = info2.strip().split("\t")
                        fw.write("%s\n" % (RE[a2[0]]))

                fr2.close()
                fw.close()
                logging.debug("Output the final result of rMATS " + a[0] + "\n")
                
fr.close()

currentTime = time.time();
runningTime = currentTime-startTime; ## in seconds
logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60));
logging.debug("Program finished")
sys.exit(0);

