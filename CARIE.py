#!/bin/python
import getopt,copy,re,os,sys,logging,time,datetime;
options, args = getopt.getopt(sys.argv[1:], 'i:o:',['input=','GTF=','output=','lib=','read=','length=','type=','anchor=','comparison=','analysis='])
input ='';
GTF ='';
output ='.'
lib ='unstrand'
read ='P'
length = 100
type = 'All'
anchor = 8
comparison = ''
analysis = 'U'
for opt, arg in options:
	if opt in ('-i','--input'):
		input = arg
	elif opt in ('--GTF'):
               	GTF = arg
	elif opt in ('-o','--output'):
		output = arg
	elif opt in ('--lib'):
		lib = arg
	elif opt in ('--read'):
                read = arg
	elif opt in ('--length'):
                length = int(arg)
	elif opt in ('--type'):
                type = arg
	elif opt in ('--anchor'):
                anchor = int(arg)
	elif opt in ('--comparison'):
               	comparison = arg
        elif opt in ('--analysis'):
                analysis = arg

if (not input or not GTF):
	print "Not enough parameters!"
	print "Program : ", sys.argv[0]
	print "          A python program to calculate the retained intron level and differential retained introns.\n"
	print "Usage :", sys.argv[0], " -i/--input: s1.bam/s1.sam[,s2.bam/s2.sam]. Mapping results for all of samples in bam/sam format. Different samples are sepreated by commas;"
	print "Usage :", sys.argv[0], " --GTF: The gtf file;"
	print "Usage :", sys.argv[0], " -o/--output: The output directory. The default is current directory;"
	print "Usage :", sys.argv[0], " --lib: The library type with choices of unstrand/first/second. The details are explained at the parameter of library-type in tophat2. The default is unstrand;"
	print "Usage :", sys.argv[0], " --read: The sequencing strategy of producing reads with choices of P/S. The default is P;"
	print "Usage :", sys.argv[0], " --length: The read length in nucleotide. The default length is 100;"
	print "Usage :", sys.argv[0], " --anchor: The anchor length in nucleotide. The program will only count reads spanning junctions with at least this anchor length on each side. The default is 8;"
	print "Usage :", sys.argv[0], " --type: The types of RI level calculation used. The choices are All/Junction/JunctionIntron/5Simple/3Simple/5Complex/3Complex. The details of each calcuation are explained in the README file. If type is All, all of the six types of calcuation will be carried out,but the rMATS will not performed. The default is All;"
	print "Usage :", sys.argv[0], " --comparison: A file providing the sample pairs needed to calculate the differential RI level.The format should be column 1(name of comparions), column 2 (sample 1 order in the input file,replicates seperated by commas), column 3 (sample 2 order in the input file,replicates seperated by commas). If absent, rMATS step will be skipped;"
	print "uasge: ", sys.argv[0], " --analysis: Type of rMATS analysis to perform. analysisType is either P or U. P is for paired analysis and U is for unpaired analysis. Default is U."
	print datetime.datetime.now()
	print "Author  : Shaofang Li"
	print "Contact : sfli001@gmail.com"
	sys.exit()

if (not os.path.exists(output)):
        os.system("mkdir %s" % output)


### setting up the logging format 
logging.basicConfig(level=logging.DEBUG,
                    format='%(asctime)s %(message)s',
                    filename=output+'/log.CAIRE.'+ str(datetime.datetime.now())+'.txt' ,
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
print path
##make directory for gtf files
gtf_path = "%s/gtf_files" % output
if (not os.path.exists(gtf_path)):
	os.system("mkdir %s" % gtf_path)
gtf = re.sub(".*/","",GTF)
cmd = "cp %s %s/%s" %(GTF,gtf_path,gtf)
os.system(cmd)
##get the path  of the python  programs
bin_path = "%s/bin" % path
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
cmd = "python %s/Clean_Intron_strand.py --gtf %s --path %s --length %s --anchor %s --strand %s" %(bin_path, gtf,gtf_path,length, anchor,lib)
os.system(cmd)
logging.debug("gtf_files\Intron_clean_" + gtf);
logging.debug("#########################################################################\n");

samples = input.split(",")
##generate the counts files
count_path = "%s/counts" % output
if (not os.path.exists(count_path)):
        os.system("mkdir %s" % count_path)
if (re.search("\.bam$",input)):
	cmd = "python %s/Intron_countall_strandBam.py --gtf %s/Intron_%s --length %s --anchor %s --bam %s -o %s/count_all.txt --lib %s --read %s" %(bin_path,gtf_path, gtf, length, anchor, ",".join(samples), count_path, lib, read)
else:
	cmd = "python %s/Intron_countall_strandSam.py --gtf %s/Intron_%s --length %s --anchor %s --sam %s -o %s/count_all.txt --lib %s --read %s" %(bin_path,gtf_path, gtf, length, anchor, ",".join(samples), count_path, lib, read)
logging.debug(cmd)
os.system(cmd)
logging.debug("finishing generating the total count file")

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
fr = open("%s/Intron_clean_%s" % (gtf_path, gtf))
for info in fr:
        a = info.strip().split("\t")
        c1= re.sub('.*clean "|\".*','',a[8])
	c5= re.sub('.*clean_5end "|\".*','',a[8])
	c3= re.sub('.*clean_3end "|\".*','',a[8])
	cc = re.sub('.*clean_simple "|\".*','',a[8])
        key = "%s_%s_%s" % (a[0],a[3],a[4])
        if(key in intron_clean):
                if(c1=="false"):
                        intron_clean[key][0]= c1
                if(c5=="false"):
                        intron_clean[key][1]= c5
                if(c3=="false"):
                        intron_clean[key][2]= c3
	else:
                intron_clean[key] = [c1,c5,c3,cc]
fr.close()

###generate the output 
output_path = "%s/result" % output
if (not os.path.exists(output_path)):
        os.system("mkdir %s" % output_path)

if(type =="Junction" or type=="All"):
	logging.debug("applying Junction method")
	fr1 =open("%s/counts/count_all.txt" % output)
	##generate the counts info for all of the intron
	fw = open("%s/counts_all_Junction.txt" %output_path, "w")
	fw.write("Intron_id\tgene_id\tstrand\tchr\tstart\tend\tannoated\tclean\tinclusion_counts\tskip_counts\tIncFormLen\tSkipFormLen\tInclusionlevel\n")
	l = len(samples)	
	for info1 in fr1:	
		a1 = info1.strip().split("\t")
		skp = [0] * l
		inc = [0] * l
		in_level=[0] *l
		sk_l = length - 2*anchor +1
                in_l = 2* (length - 2*anchor +1 )
		for i in range(0,l):
			inc[i] = str(int(a1[i*6+6]) + int(a1[i*6+8]))
			skp[i] = a1[i*6+10]
			if((inc[i] + skp[i]) > 0) & ((float(inc[i])/in_l + float(skp[i])/sk_l) > 0):	
				in_level[i]= str((float(inc[i])/in_l)/(float(inc[i])/in_l + float(skp[i])/sk_l))
			else:
				in_level[i] ='NA'
		fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (a1[0],a1[1],a1[2],a1[3], a1[4],a1[5],intron_anno[a1[0]],",".join(intron_clean[a1[0]]), ",".join(inc), ",".join(skp), in_l,sk_l,",".join(in_level)))
	fw.close()
	fr1.close()
	logging.debug("done generating the count file")	

if(type =="5Simple" or type =="All"):
        logging.debug("applying 5Simple method")
	fr1 =open("%s/counts/count_all.txt" % output)
        ##generate the counts info for all of the intron
        fw = open("%s/counts_all_5Simple.txt" %output_path, "w")
        fw.write("Intron_id\tgene_id\tstrand\tchr\tstart\tend\tannoated\tclean\tinclusion_counts\tskip_counts\tIncFormLen\tSkipFormLen\tInclusionlevel\n")
	l = len(samples)
        for info1 in fr1:
                a1 = info1.strip().split("\t")
                skp = [0] * l
                inc = [0] * l
		in_level=[0] *l
                sk_l = length - 2*anchor +1
                in_l = length - 2*anchor +1
		for i in range(0,l):
			inc[i] = a1[i*6+6]
			skp[i] = a1[i*6+10]
			if((inc[i] + skp[i]) > 0) & ((float(inc[i])/in_l + float(skp[i])/sk_l) > 0):
                                in_level[i]= str((float(inc[i])/in_l)/(float(inc[i])/in_l + float(skp[i])/sk_l))
                        else:
                                in_level[i] ='NA'
		fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (a1[0],a1[1],a1[2],a1[3], a1[4],a1[5],intron_anno[a1[0]],",".join(intron_clean[a1[0]]), ",".join(inc), ",".join(skp), in_l,sk_l,",".join(in_level)))
	fw.close()
        fr1.close()
	logging.debug("done generating the count file") 

if(type =="3Simple"or type=="All"):
	logging.debug("applying 3Simple method")
	fr1 =open("%s/counts/count_all.txt" % output)
        ##generate the counts info for all of the intron
        fw = open("%s/counts_all_3Simple.txt" %output_path, "w")
        fw.write("Intron_id\tgene_id\tstrand\tchr\tstart\tend\tannoated\tclean\tinclusion_counts\tskip_counts\tIncFormLen\tSkipFormLen\tInclusionlevel\n")
        l = len(samples)
        for info1 in fr1:
                a1 = info1.strip().split("\t")
                skp = [0] * l 
                inc = [0] * l 
                in_level=[0] *l
                sk_l = length - 2*anchor +1
                in_l = length - 2*anchor +1
                for i in range(0,l):
                        inc[i] = a1[i*6+8]
                        skp[i] = a1[i*6+10]
                  	if((inc[i] + skp[i]) > 0) & ((float(inc[i])/in_l + float(skp[i])/sk_l) > 0):
                                in_level[i]= str((float(inc[i])/in_l)/(float(inc[i])/in_l + float(skp[i])/sk_l))
                        else:
                                in_level[i] ='NA'
                fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (a1[0],a1[1],a1[2],a1[3], a1[4],a1[5],intron_anno[a1[0]],",".join(intron_clean[a1[0]]), ",".join(inc), ",".join(skp), in_l,sk_l,",".join(in_level)))
        fw.close()
        fr1.close()
        logging.debug("done generating the count file")


if(type =="5Complex" or type=="All"):
        logging.debug("applying 5Complex method")
	fr1 =open("%s/counts/count_all.txt" % output)
        ##generate the counts info for all of the intron
        fw = open("%s/counts_all_5Complex.txt" %output_path, "w")
        fw.write("Intron_id\tgene_id\tstrand\tchr\tstart\tend\tannoated\tclean\tinclusion_counts\tskip_counts\tIncFormLen\tSkipFormLen\tInclusionlevel\n")
        l = len(samples)
        for info1 in fr1:
                a1 = info1.strip().split("\t")
                skp = [0] * l 
                inc = [0] * l 
                in_level=[0] *l
                sk_l = length - 2*anchor +1
                in_l = length - 2*anchor +1
                for i in range(0,l):
                        inc[i] = a1[i*6+6]
                        skp[i] = a1[i*6+7]
                	if((inc[i] + skp[i]) > 0) & ((float(inc[i])/in_l + float(skp[i])/sk_l) > 0):
                                in_level[i]= str((float(inc[i])/in_l)/(float(inc[i])/in_l + float(skp[i])/sk_l))
                        else:
                                in_level[i] ='NA'
		fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (a1[0],a1[1],a1[2],a1[3], a1[4],a1[5],intron_anno[a1[0]],",".join(intron_clean[a1[0]]), ",".join(inc), ",".join(skp), in_l,sk_l,",".join(in_level)))
        fw.close()
        fr1.close()
        logging.debug("done generating the count file")

if(type =="3Complex" or type =="All"):
	logging.debug("applying 3Complex method")
	fr1 =open("%s/counts/count_all.txt" % output)
        ##generate the counts info for all of the intron
        fw = open("%s/counts_all_3Complex.txt" %output_path, "w")
        fw.write("Intron_id\tgene_id\tstrand\tchr\tstart\tend\tannoated\tclean\tinclusion_counts\tskip_counts\tIncFormLen\tSkipFormLen\tInclusionlevel\n")
        l = len(samples)
        for info1 in fr1:
                a1 = info1.strip().split("\t")
                skp = [0] * l
                inc = [0] * l 
                in_level=[0] *l
                sk_l = length - 2*anchor +1
                in_l = length - 2*anchor +1
                for i in range(0,l):
                        inc[i] = a1[i*6+8]
                        skp[i] = a1[i*6+9]
			if((inc[i] + skp[i]) > 0) & ((float(inc[i])/in_l + float(skp[i])/sk_l) > 0):
                                in_level[i]= str((float(inc[i])/in_l)/(float(inc[i])/in_l + float(skp[i])/sk_l))
                        else:
                                in_level[i] ='NA'        
                fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (a1[0],a1[1],a1[2],a1[3], a1[4],a1[5],intron_anno[a1[0]],",".join(intron_clean[a1[0]]), ",".join(inc), ",".join(skp), in_l,sk_l,",".join(in_level)))
        fw.close()
        fr1.close()
        logging.debug("done generating the count file")

if(type =="JunctionIntron" or type =="All"):
	logging.debug("applying JunctionIntron method")
	fr1 =open("%s/counts/count_all.txt" % output)
        ##generate the counts info for all of the intron
        fw = open("%s/counts_all_JunctionIntron.txt" %output_path, "w")
        l = len(samples)
	fw.write("Intron_id\tgene_id\tstrand\tchr\tstart\tend\tannoated\tclean\tinclusion_counts\tskip_counts\tIncFormLen\tSkipFormLen\tInclusionlevel\n")
	for info1 in fr1:
                a1 = info1.strip().split("\t")
                skp = [0] * l
                inc = [0] * l
                in_level=[0] *l
                sk_l = length - 2*anchor +1
                in_l = length - 2*anchor +2 + int(a1[5])-int(a1[4])
                for i in range(0,l):
                        inc[i] =str(int(a1[i*6+6]) + int(a1[i*6+8])+ int(a1[i*6+11]))
                        skp[i] = a1[i*6+10]
                	if((inc[i] + skp[i]) > 0) & ((float(inc[i])/in_l + float(skp[i])/sk_l) > 0):
                                in_level[i]= str((float(inc[i])/in_l)/(float(inc[i])/in_l + float(skp[i])/sk_l))
                        else:
                                in_level[i] ='NA'
		fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (a1[0],a1[1],a1[2],a1[3], a1[4],a1[5],intron_anno[a1[0]],",".join(intron_clean[a1[0]]), ",".join(inc), ",".join(skp), in_l,sk_l,",".join(in_level)))

	
        fw.close()
        fr1.close()
	logging.debug("done generating the count file")


if (not os.path.exists(comparison) or not comparison):
	logging.debug("The comparion file was not provided")
	currentTime = time.time();
        runningTime = currentTime-startTime; ## in seconds
        logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60));
	logging.debug("Program finished")
	sys.exit()
if(type=="All"):
	logging.debug("The type of count file was not specified")
	currentTime = time.time();
	runningTime = currentTime-startTime; ## in seconds
	logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60));
	logging.debug("Program finished")
	sys.exit()
logging.debug("running rMATS ")
fr =open(comparison)
for info in fr:
        a = info.strip().split()
	ss1= a[1].split(",")
        ss2= a[2].split(",")
	fw = open("%s/counts/rMATS_%s_%s.txt" % (output,a[0],type),"w")	
	fw.write("ID\tIJC_SAMPLE_1\tSJC_SAMPLE_1\tIJC_SAMPLE_2\tSJC_SAMPLE_2\tIncFormLen\tSkipFormLen\n")
	fr1 = open("%s/counts_all_%s.txt" % (output_path,type))
	info1 = fr1.readline()
	for info1 in fr1:
		a1 = info1.split("\t")
		skp1 = [0] * len(ss1)
		inc1 = [0] * len(ss1)
		skp2 = [0] * len(ss2)
		inc2 = [0] * len(ss2)
		in_level = a1[8].split(",")
		sk_level = a1[9].split(",")
		for i in range(0,len(ss1)):
			inc1[i] = in_level[int(ss1[i])-1]
			skp1[i] = sk_level[int(ss1[i])-1]
		for i in range(0,len(ss2)):
                        inc2[i] = in_level[int(ss2[i])-1]
                        skp2[i] = sk_level[int(ss2[i])-1]	
		
		fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (a1[0],",".join(inc1),",".join(skp1), ",".join(inc2),",".join(skp2),a1[10],a1[11]))
	fr1.close()
	fw.close()
	cmd ="%s/MATS/rMATS.sh -d %s/counts/rMATS_%s_%s.txt -o %s/counts/rMATS_%s_%s -p 1  -t %s -c 0.0001" %(bin_path,output,a[0],type,output,a[0],type,analysis)
        os.system(cmd)
	logging.debug("done running the rMATS for" + a[0])
	logging.debug("output the final result of rMATS" + a[0])
	fr1 = open("%s/counts_all_%s.txt" % (output_path,type))
	fr2 = open("%s/counts/rMATS_%s_%s/rMATS_Result.txt" %(output,a[0],type))
	fw = open("%s/result/rMATS_Result_%s_%s.txt" %(output, a[0],type),"w")
	fw.write("Intron_id\tgene_id\tstrand\tchr\tstart\tend\tannoated\tclean\tinclusions_counts_SAMPLE1\tskip_counts_SAMPLE1\tinclusions_counts_SAMPLE2\tskip_counts_SAMPLE2\tInclusion_length\tSkipping_length\tPValue\tFDR\tIncLevel_SAMPLE1\tIncLevel_SAMPLE2\tIncLevelDifference\n")
	info1 = fr1.readline()
        info2 = fr2.readline()
	info1 = fr1.readline()
        info2 = fr2.readline()
	while(info1):
                a1 = info1.strip().split("\t")
                a2 = info2.strip().split("\t")
                fw.write("%s\t%s\n" % ("\t".join(a1[0:8]), "\t".join(a2[1:])))
                info1 = fr1.readline()
                info2 = fr2.readline()
        fr1.close()
        fr2.close()
        fw.close()
fr.close()

currentTime = time.time();
runningTime = currentTime-startTime; ## in seconds
logging.debug("Program ran %.2d:%.2d:%.2d" % (runningTime/3600, (runningTime%3600)/60, runningTime%60));
logging.debug("program finished")
sys.exit(0);

