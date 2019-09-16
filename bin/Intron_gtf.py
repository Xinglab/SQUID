#!/usr/bin/python
import getopt,re,os,sys,logging,time,datetime,copy;
options, args = getopt.getopt(sys.argv[1:],'', ['gtf=','path='])
gtf = ''
path=''
for opt, arg in options:
        if opt in('--gtf'):
                gtf = arg
        elif opt in('--path'):
                path = arg
if (not gtf or not path):
        print "not enough parameters"
        print "Program : ", sys.argv[0]
        print "a python program to get the intron from gtf file with flanking exons info available"
        print "usage :", sys.argv[0], " --gtf: gtf file"
        print "usage :", sys.argv[0], " --path: the directory to the input and output gtf file"
        print datetime.datetime.now()
        print "Author  : Shaofang Li"
        print "Contact : sfli001@gmail.com"
        sys.exit()
fr = open ("%s/%s" % (path,gtf))
exon_s=dict()
exon_e=dict()
info1 = fr.readline()
while (info1[0]=="#"):
	info1= fr.readline()
a1 = info1.strip().split("\t")
while (a1[2]!="exon"):
        info1 = fr.readline()
	a1 = info1.strip().split("\t")
while info1:
	a1 = info1.strip().split("\t")
        if (a1[2]!= "exon"):
                info1=fr.readline()
                continue
	transcript_id = re.sub('.*transcript_id "|\".*','',a1[8])
#	gene_id = re.sub('.*gene_id "|\".*','',a1[8])
#	gene_name = re.sub('.*gene_name "|\".*','',a1[8])
	exon_s.setdefault(transcript_id,[]).append(int(a1[3])) 
	exon_e.setdefault(transcript_id,[]).append(int(a1[4]))
	info1 = fr.readline()
fr.close()

#print "exon info"
fr = open ("%s/%s" % (path,gtf))
fw = open(("%s/Intron_%s" %(path, gtf)),"w")
info1 = fr.readline()
while (info1[0]=="#"):
        info1= fr.readline()
a1 = info1.strip().split("\t")
while (a1[2]!="exon"):
        info1 = fr.readline()
        a1 = info1.strip().split("\t")
trans = dict()
while info1:
        a1 = info1.strip().split("\t")
	if(a1[2] !="exon"):
		info1 = fr.readline()	
		continue
	transcript_id = re.sub('.*transcript_id "|\".*','',a1[8])
	gene_id = re.sub('.*gene_id "|\".*','',a1[8])
        gene_name = re.sub('.*gene_name "|\".*','',a1[8])
	
        if(transcript_id in trans):
		trans[transcript_id]="false"
	elif(len(exon_s[transcript_id])>1):
		exon_s[transcript_id].sort()
		exon_e[transcript_id].sort()
		for i in range(1,len(exon_s[transcript_id])):
			if(a1[6]=="+"):
				fw.write("%s\t%s\tintron\t%s\t%s\t.\t%s\t.\tgene_id \"%s\"; gene_name \"%s\"; transcript_id \"%s\"; intron_number \"%s\"; total_intron_number \"%s\"; exon_1 \"%s_%s\"; exon_2 \"%s_%s\";\n"%(a1[0],a1[1],str(int(exon_e[transcript_id][i-1])+1),str(int(exon_s[transcript_id][i])-1),a1[6],gene_id,gene_name,transcript_id,i,len(exon_s[transcript_id])-1,exon_s[transcript_id][i-1],exon_e[transcript_id][i-1],exon_s[transcript_id][i],exon_e[transcript_id][i]))
			else:
				fw.write("%s\t%s\tintron\t%s\t%s\t.\t%s\t.\tgene_id \"%s\"; gene_name \"%s\"; transcript_id \"%s\"; intron_number \"%s\"; total_intron_number \"%s\"; exon_1 \"%s_%s\"; exon_2 \"%s_%s\";\n"%(a1[0],a1[1],str(int(exon_e[transcript_id][i-1])+1),str(int(exon_s[transcript_id][i])-1),a1[6],gene_id,gene_name,transcript_id,len(exon_s[transcript_id])-i,len(exon_s[transcript_id])-1,exon_s[transcript_id][i-1],exon_e[transcript_id][i-1],exon_s[transcript_id][i],exon_e[transcript_id][i]))
		trans[transcript_id]="true"
	info1 = fr.readline()
	
fr.close()
fw.close()
