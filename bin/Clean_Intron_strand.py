#!/bin/python
import copy, getopt,re,os,sys,logging,time,datetime;

options, args = getopt.getopt(sys.argv[1:],'', ['gtf=','path=','length=','anchor=','strand='])
gtf = ''
path=''
strand = 'unstrand'
length =100
anchor =8
for opt, arg in options:
        if opt in('--gtf'):
                gtf = arg
        elif opt in('--path'):
                path = arg
	elif opt in('--length'):
                length = int(arg)
	elif opt in('--anchor'):
                anchor = int(arg)
	elif opt in('--strand'):
                strand = arg
if (not gtf or not path):
        print "Not enough parameters!"
        print "Program : ", sys.argv[0]
        print "          A python program to get the clean intron  from given intron gtf file."
        print "Usage :", sys.argv[0], " --gtf: The gtf file;"
        print "Usage :", sys.argv[0], " --path: The directory to the input and output gtf file;"
	print "Usage :", sys.argv[0], " --length: The read length;"
	print "Usage :", sys.argv[0], " --anchor: The anchor length;"
        print "Usage :", sys.argv[0], " --strand: The library type."
	print datetime.datetime.now()
        print "Author  : Shaofang Li"
        print "Contact : sfli001@gmail.com"
        sys.exit()

fr = open ("%s/Intron_%s" % (path,gtf))
fw = open(("%s/Intron_clean_%s" %(path,gtf)),"w")



pos = dict()
bin = 1000
intron = dict()
intron5= dict()
intron3 = dict()
for info1 in fr:
	a1 = info1.strip().split("\t")
	key = "%s_%s_%s" % (a1[0],a1[3],a1[4])
	intron[key]=["true","true","true","true",a1[0],int(a1[3]),int(a1[4])]
	intron5.setdefault((a1[0],int(a1[3])),[]).append(int(a1[4]))	
        intron3.setdefault((a1[0],int(a1[4])),[]).append(int(a1[3]))
	index_s =  int(a1[3])/ bin
        index_e =  int(a1[4])/ bin
	if(strand=="unstrand"):
		for i in  range(index_s , index_e+1):
                        pos.setdefault((a1[0],i),[]).append(key)
	else:
        	for i in  range(index_s , index_e+1):
                	pos.setdefault((a1[0],a1[6],i),[]).append(key)
		
fr.close()
for k in intron5:
	intron5[k] = list(set(intron5[k]))
for k in intron3:
        intron3[k] = list(set(intron3[k]))

for k in intron:
	for ee in intron5[intron[k][4],intron[k][5]]:
		if(ee!=intron[k][6]):
			intron[k][3]="false"
	for ss in intron3[intron[k][4],intron[k][6]]:
		if(ss!=intron[k][5]):
			intron[k][3]="false"

fr2 = open ("%s/Exon_%s" % (path,gtf))
for info2 in fr2:
        a2 = info2.strip().split("\t")
	gene_id = re.sub('.*gene_id "|\".*','',a2[8])
        index_s =  int(a2[3]) / bin
        index_e =  int(a2[4]) / bin
	if(strand=="unstrand"):
                for i in  range(index_s , index_e+1):
                        if (a2[0], i) in pos:
                                        for j in pos[a2[0], i]: 
                                                if(int(a2[3]) <= intron[j][5]):
                                                        if(int(a2[4]) >= intron[j][5]):
                                                                intron[j][0]="false"
                                                                intron[j][1]="false"
                                                if( (int(a2[3]) > intron[j][5]) and  (int(a2[3]) <= intron[j][6])):
                                                        intron[j][0]="false"
                                                if( (int(a2[3]) > intron[j][5]) and (int(a2[3]) <= min(intron[j][6],intron[j][5]+ length -anchor +1))):
                                                        intron[j][1]="false"
                                                        if(intron[j][0]=="true"):
                                                                print info2, intron[j]
                                                if(int(a2[3]) <= max(intron[j][5], intron[j][6]-length+ anchor+1)):
                                                        if(int(a2[4]) >= max(intron[j][5], intron[j][6]-length+ anchor+1)):
                                                                intron[j][2]="false"
                                                if( (int(a2[3]) > max(intron[j][5],intron[j][6]-length+anchor +1)) and  (int(a2[3]) <= intron[j][6])):                                                        intron[j][2]="false"
                                                                        
                                                if(intron[j][0] =="true" and intron[j][2]=="false"):
                                                        print info2, intron[j]
	else:	
		for i in  range(index_s , index_e+1):
			if (a2[0], a2[6], i) in pos:
					for j in pos[a2[0],a2[6], i]:
						if(int(a2[3]) <= intron[j][5]):
							if(int(a2[4]) >= intron[j][5]):
								intron[j][0]="false"
								intron[j][1]="false"
						if( (int(a2[3]) > intron[j][5]) and  (int(a2[3]) <= intron[j][6])):
							intron[j][0]="false"
						if( (int(a2[3]) > intron[j][5]) and (int(a2[3]) <= min(intron[j][6],intron[j][5]+ length -anchor +1))):
							intron[j][1]="false"
							if(intron[j][0]=="true"):
								print info2, intron[j]
						if(int(a2[3]) <= max(intron[j][5], intron[j][6]-length+ anchor+1)):
							if(int(a2[4]) >= max(intron[j][5], intron[j][6]-length+ anchor+1)):
								intron[j][2]="false"
						if( (int(a2[3]) > max(intron[j][5],intron[j][6]-length+anchor +1)) and  (int(a2[3]) <= intron[j][6])):
							intron[j][2]="false"	
									
						if(intron[j][0] =="true" and intron[j][2]=="false"):
							print info2, intron[j]		
						
					

fr2.close()

fr = open ("%s/Intron_%s" % (path,gtf))
for info1 in fr:
	a1 = info1.strip().split("\t")
	key = "%s_%s_%s" % (a1[0],a1[3],a1[4])
	if(a1[6]=="+"):
		fw.write('%s clean "%s"; clean_5end "%s"; clean_3end "%s"; clean_simple "%s";\n' %(info1.strip(),intron[key][0],intron[key][1],intron[key][2], intron[key][3]))
	if(a1[6]=="-"):
		fw.write('%s clean "%s"; clean_5end "%s"; clean_3end "%s"; clean_simple "%s";\n' %(info1.strip(),intron[key][0],intron[key][2],intron[key][1],intron[key][3]))
		
fr.close()
fw.close()

