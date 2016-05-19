#!/bin/python
import getopt,copy,re,os,sys,logging,time,datetime;
import pysam,os.path
options, args = getopt.getopt(sys.argv[1:], 'i:o:',['input=','gtf=','exon=','lim=','lib=','output='])
input='';
gtf=''
exon ='';
lim=5;
lib = 'unstrand'
output='';
for opt, arg in options:
	if opt in ('-i','--input'):
		input = arg
	elif opt in ('--gtf'):
                gtf = arg
	elif opt in ('--exon'):
		exon = arg
	elif opt in ('--lim'):
                lim = int(arg)
	elif opt in ('--lib'):
                lib = arg
	elif opt in ('-o','--output'):
                output = arg
if (not input or not exon or not output or not gtf):
	print "Not enough parameters!"
	print "Program : ", sys.argv[0]
	print "          A python program to get the adjusted the intron based on the junction count and exon gtf file"
	print "Usage :", sys.argv[0], " --input: the junction input file;"
	print "Usage :", sys.argv[0], " --gtf: the intron gtf file with annotation of clean intron or not"
	print "Usage :", sys.argv[0], " --exon: exon gtf file "
	print "Usage :", sys.argv[0], " --lim: the number of the read threshhold to be considered"
	print "Usage :", sys.argv[0], " --lib: the library type"
	print "Usage :", sys.argv[0], ' --output: first column intron ID, second, chr, third, adjusted start, fourth adjusted end'
	print datetime.datetime.now()
	print "Author  : Shaofang Li"
	print "Contact : sfli001@gmail.com"
	sys.exit()
##library info was not directly used in this program, the strand info was utilized in the junction and clean intron files
bin = 1000
intron = dict()
intronJ = dict()
intronE = dict()
pos_in = dict()
fr = open(gtf)
for info in fr:
	a = info.strip().split("\t")
	C1= re.sub('.*clean "|\".*','',a[8])
	if(C1 =="false"):
		key ="%s_%s_%s" % (a[0],a[3],a[4])
		if( not intron.has_key(key)):
			intron[key]= [a[0],int(a[3]),int(a[4])]	
			intronE[key] =[]
			index_s =  int(a[3]) / bin
			index_e =  int(a[4]) / bin
			for i in  range(index_s , index_e+1):
				pos_in.setdefault((a[0],i),[]).append(key)
				
fr.close()
##generate the dictionary with junction read
fr = open(input)	
for info in fr:
	a = info.strip().split("\t")
	if(int(a[1]) < lim):
		continue
	jun_p = a[0].split("_")
	for j in range(2, len(a)):
		intronJ.setdefault(a[j],[]).append((jun_p[0],int(jun_p[1]),int(jun_p[2])))
		if( not intron.has_key(a[j])):
			in_p = a[j].split("_")
			index_s =  int(in_p[1]) / bin
			index_e =  int(in_p[2]) / bin
			intron[a[j]]= [in_p[0], int(in_p[1]),int(in_p[2])]
			for i in  range(index_s , index_e+1):
				pos_in.setdefault((a[0],i),[]).append(key)
fr.close()
##generate the dictionary with overlapping exons
fr= open(exon)
exon = dict()
for info in fr:
	a = info.strip().split("\t")
	key = "%s_%s_%s" % (a[0],a[3],a[4])
	if(exon.has_key(key)):
		continue
	exon[key]="true"
	index_s =  int(a[3]) / bin
	index_e =  int(a[4]) / bin
	for i in  range(index_s , index_e+1):
		if(not pos_in.has_key((a[0],i))):
			continue
		for ii in pos_in[a[0],i]:
			in_p = ii.split("_")
			if(int(in_p[1]) > int(a[4])):
				continue
			elif( int(in_p[2]) < int(a[3])):
				continue
			else:
				intronE.setdefault(ii,[]).append((a[0],int(a[3]), int(a[4])))
	
fr.close()
#concatenate the exon and junction region
SS = 0
for key in intron:
	in_p = key.split("_")
	if(key in intronE ):
		SS +=1 	
		intronE[key] = list(set(intronE[key]))
	
		if(len(intronE[key]) !=1):
			exonS = []
			exond = dict()
			for ee in intronE[key]:
				exonS.append(ee[1])
				if(exond.has_key(ee[1]) and int(ee[2]) > exond[ee[1]]):
					exond[ee[1]]= int(ee[2])
				if(not exond.has_key(ee[1])):
					exond[ee[1]]= int(ee[2])
			exonS = list(set(exonS))
			if(len(exonS) ==1):
				intronE[key] =[]
				intronE[key].append((in_p[0],exonS[0],exond[exonS[0]]))
			else:
				exonS.sort()
				exonS1 = []
				exond1 = dict()	
				s1 = exonS[0]
				e1 = exond[exonS[0]]
				for i in range(1,len(exonS)):
					s2 = exonS[i]
					st = min(s1,s2)
					e2 = exond[exonS[i]]
					end = max(e1, e2)
					if((end -st -1) <= (e1-s1 + e2 -s2)):
						s1 = st
						e1 = end
						
					else:
						exonS1.append(s1)
						exond1[s1] = e1
						s1 = s2
						e1 = e2
					if(i == len(exonS)-1):
						exonS1.append(s1)
						exond1[s1]=e1
				intronE[key]=[]
				for st in exonS1:
					intronE[key].append((in_p[0],st,exond1[st]))
	#		print "intronE adjust", intronE[key]
	if(key in intronJ and len(intronJ[key])==1):
		JUN = intronJ[key][0]
		intronJ[key]= [(in_p[0],max(JUN[1],int(in_p[1])),min(int(in_p[2]),JUN[2]))]
	if(key in intronJ and len(intronJ[key]) > 1):
		junS = []
		jund = dict()
		for jun in intronJ[key]:
			junS.append(jun[1])
			if(jund.has_key(jun[1]) and jun[2] <  jund[jun[1]]):
				jund[jun[1]]= jun[2]
			if(not jund.has_key(jun[1])):
				jund[jun[1]]= jun[2]
		junS = list(set(junS))
		if(len(junS) ==1):
			intronJ[key] =[]
			intronJ[key].append((in_p[0],max(junS[0],int(in_p[1])),min(int(in_p[2]),jund[junS[0]])))
		else:	
			junS.sort()
			junS1 = []
			jund1 = dict()
			s1 = junS[0]
			e1 = jund[junS[0]]
			for i in range(1,len(junS)):
				s2 = junS[i]
				st = max(s1,s2)
				e2 = jund[junS[i]]
				end = min(e1, e2)
				if(end > st):
					s1 = st
					e1 = end
				else:
					junS1.append(s1)
					jund1[s1] = e1
					s1 = s2
					e1 = e2
				if(i == len(junS)-1):
					junS1.append(s1)
					jund1[s1]=e1
			intronJ[key]=[]
			for st in junS1:
				intronJ[key].append((in_p[0] ,max(int(in_p[1]),st),min(int(in_p[2]),jund1[st])))
##get the optimal region using EXON and junction reads
IntronA = dict()
for key in intron:
	in_p = key.split("_")
	if((key in intronE) and (key not in intronJ)):
		if(int(in_p[1]) < intronE[key][0][1]-1):
			IntronA.setdefault(key,[]).append("%s_%s_%s" % (in_p[0], in_p[1],int(intronE[key][0][1])-1))
		for i in range(1, len(intronE[key])):	
			E = int(intronE[key][i][1]) -1
			S = int(intronE[key][i-1][2])+ 1
			IntronA.setdefault(key,[]).append("%s_%s_%s" % (in_p[0], S,E))
		if(int(in_p[2]) > int(intronE[key][len(intronE[key])-1][2])+1):
			IntronA.setdefault(key,[]).append("%s_%s_%s" % (in_p[0], int(intronE[key][len(intronE[key])-1][2])+1,in_p[2]))
		if(not IntronA.has_key(key)):
			IntronA.setdefault(key,[]).append("%s_%s_%s" % (in_p[0], in_p[1], in_p[1]))
	if((key not in intronE) and (key in intronJ)):
		for i in range(0, len(intronJ[key])):
			IntronA.setdefault(key,[]).append("%s_%s_%s"% (intronJ[key][i][0], intronJ[key][i][1],intronJ[key][i][2]))
	if((key in intronE) and (key in intronJ)):
		JUN1 =[]
		EXON = intronE[key]
		j = 0
		for i in range(0, len(intronJ[key])):
			if(int(EXON[j][1] )> int(intronJ[key][i][2])):
				IntronA.setdefault(key,[]).append("%s_%s_%s"% (intronJ[key][i][0], intronJ[key][i][1],intronJ[key][i][2]))
				continue
			while(j < (len(intronE[key])-1) and int(EXON[j][2])< int(intronJ[key][i][1])):
				j+=1
			if(EXON[j][2] < intronJ[key][i][1]):
				IntronA.setdefault(key,[]).append("%s_%s_%s"% (intronJ[key][i][0], intronJ[key][i][1],intronJ[key][i][2]))
				continue
			
			if(EXON[j][2] > intronJ[key][i][2] and EXON[j][1] > intronJ[key][i][1]):
				IntronA.setdefault(key,[]).append("%s_%s_%s"% (intronJ[key][i][0], intronJ[key][i][1],EXON[j][1]-1))
			if(EXON[j][2] < intronJ[key][i][2] and EXON[j][1] >  intronJ[key][i][1]):
				IntronA.setdefault(key,[]).append("%s_%s_%s"% (intronJ[key][i][0], intronJ[key][i][1],EXON[j][1]-1))
			if(EXON[j][2] < intronJ[key][i][2]):
				if(j ==len(EXON)-1):
					IntronA.setdefault(key,[]).append("%s_%s_%s"% (intronJ[key][i][0], EXON[j][2]+1,intronJ[key][i][2]))
					continue
				while (j < (len(intronE[key])-1) and EXON[j+1][1] -1  < intronJ[key][i][2]):
					IntronA.setdefault(key,[]).append("%s_%s_%s"% (intronJ[key][i][0], EXON[j][2]+1,EXON[j+1][1]-1))
					j+=1
		if(not IntronA.has_key(key)):
			IntronA.setdefault(key,[]).append("%s_%s_%s" % (in_p[0], in_p[1], in_p[1]))	
#print len(intron)
#print len(intronJ)
# print len(intronE)
#	print len(IntronA)
fw = open(output,"w")
for key in IntronA:
	fw.write("%s\t%s\n" % (key, "\t".join(IntronA[key])))
fw.close() 
