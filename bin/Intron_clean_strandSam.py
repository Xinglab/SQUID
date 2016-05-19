#!/bin/python
import getopt,copy,re,os,sys,logging,time,datetime;
import pysam,os.path
options, args = getopt.getopt(sys.argv[1:], 'i:o:',['input=','countfile=','anchor=','lib=','read=','length=','sam=','output=',])
input='';
read='';
sam='';
lib='';
length =0
anchor =0
output='';
countfile='';
for opt, arg in options:
	if opt in ('-o','--output'):
		output = arg
	elif opt in ('-i','--input'):
		input = arg
	elif opt in ('--countfile'):
                countfile = arg
	elif opt in ('--sam'):
                sam = arg
	elif opt in ('--lib'):
                lib = arg
	elif opt in ('--read'):
                read = arg
	elif opt in ('--length'):
                length= int(arg)
	elif opt in ('--anchor'):
                anchor= int(arg)
if (not input or not read or not sam or not output or not length or not anchor):
	print "Not enough parameters!"
	print "Program : ", sys.argv[0]
	print "          A python program to count the reads for retained intron events for varities of junction from a series of sam file."
	print "Usage :", sys.argv[0], " --input: the intron input file;"
	print "Usage :", sys.argv[0], " --countfile: the count file generated at the first step;"
	print "Usage :", sys.argv[0], " --length:the length of reads;"
	print "Usage :", sys.argv[0], " --anchor:the anchor length of the read;"
	print "Usage :", sys.argv[0], " --sam: the sam file,multiple sam file seperated by commas;"
	print "Usage :", sys.argv[0], " --lib: the library type;"
	print "Usage :", sys.argv[0], " --read: The sequencing strategy of producing reads with choices of P/S;"
	print "Usage :", sys.argv[0], ' --output: intron_id, gene_id,strand,chr,start,end,5SS inclusion counts,5SS skipping counts,3SS includion counts,3SS skipping counts,skipping counts,intron counts.'
	print datetime.datetime.now()
	print "Author  : Shaofang Li"
	print "Contact : sfli001@gmail.com"
	sys.exit()
if(lib=="unstrand"):
	intron = dict()
	fr = open(countfile)
	for info in fr:
		a = info.strip().split("\t")
		intron[a[0]] = [info]
	fr.close() 
	count = dict()
	#the dict to store the intron count,first is the inclusion at left side, second is the inclusion at right side, third is the intron counts
	ppL = dict()
        #the dict to store left edge  position of intron
        ppR = dict()
        #the dict to store right  edge  position of intron
	pos = dict()
	# the dictionary to store the introns in a certain windows
	bin = 1000
	# set the bin size
	samfile = sam.split(",")
	num = len(samfile)
	fr1 = open(input)
	for info1 in fr1:
		a1 = info1.strip().split("\t")
		for key in a1[1:]:
			intron[a1[0]].append(key)
			if(key in count):
				continue
			if(key in intron):
				continue
			count[key]= [0]*num*3
			count[key].append("true")
			in_p = key.split("_")
			ppL.setdefault((in_p[0],int(in_p[1])),[]).append(key)
			ppR.setdefault((in_p[0],int(in_p[2])),[]).append(key)
			index_s =  int(in_p[1]) / bin
			index_e =  int(in_p[2]) / bin
			for i in  range(index_s , index_e+1):
				pos.setdefault((in_p[0],i),[]).append(key)
			
	fr1.close()
	samfile = sam.split(",")
	nn = 0
	for file in samfile:
		fr2 = open(file)
		for info2 in fr2:
			if (info2[0]== "@"):
                                continue
			a2 = info2.strip().split("\t")
			if(read =="P"):
                                if((int(a2[1])/2 )%2 ==0):
                                        continue
                                if(re.search("NH\:i\:1[0-9]",info2)):
                                        continue
                                if(re.search("NH\:i\:[^1]",info2)):
                                        continue
                                if(re.search("D|H|S",a2[5])):
                                        continue
                                if(re.search("I",a2[5])):
                                        continue
                        if(read =="S"):
                                if(re.search("NH\:i\:1[0-9]",info2)):
                                        continue
                                if(re.search("NH\:i\:[^1]",info2)):
                                        continue
                                if(re.search("D|H|S",a2[5])):
                                        continue
                                if(re.search("I",a2[5])):
                                        continue

                        aa1 = a2[5].split("M")
                        aa2 = a2[5].split("N")
                        l = len(aa2)
			if(len(aa1)== 2):
				for p in range(int(a2[3]) + anchor, int(a2[3]) + int(aa1[0]) - anchor +1):
                                        if(a2[2],p) in ppL:
                                                for id in ppL[a2[2],p]:
                                                        #inclusion count at left side
                                                        count[id][nn*3] +=1
                                                        continue

				for p in range(int(a2[3]) + anchor-1, int(a2[3]) + int(aa1[0]) - anchor):
                                        if(a2[2],p) in ppR:
                                                for id in ppR[a2[2],p]:
                                                        #inclusion count at right side
                                                        count[id][nn*3+1] +=1 
                                                        continue        
                                index= int(a2[3])/bin
                                if (a2[2], index) in pos:
                                        for key in pos[a2[2],index]:
                                                in_p = key.split("_")
                                                if ((int(in_p[1]) <= int(a2[3]) )& (int(in_p[2])>= (int(a2[3]) + length -1))):
                                                        count[key][nn*3+2]+=1
			start = int(a2[3])
                        for i in range(0,l-1):
                                if(re.search("\D",aa2[i].split("M")[0])):
                                        continue
                                if(re.search("\D",aa2[i].split("M")[1])):
                                        continue
                                if(re.search("\D",aa2[i+1].split("M")[0])):
                                        continue
                                n1=int(aa2[i].split("M")[0])
                                n2=int(aa2[i].split("M")[1])
                                n3=int(aa2[i+1].split("M")[0])
				if( n1 -anchor+1 >  anchor):
                                        ss = start
                                        for p in range( ss + anchor, ss + n1 -anchor+1 ):
                                                if(a2[2],p) in ppL:
                                                        for id in ppL[a2[2],p]:
                                                                #included count at left side
                                                                count[id][nn*3] +=1
                                        for p in range( ss + anchor-1, ss + n1 -anchor ):
                                                if(a2[2],p) in ppR:
                                                        for id in ppR[a2[2],p]:
                                                                #included count at right side
                                                                count[id][nn*3+1] +=1
                                if( n3 -anchor+1 >  anchor):
                                        ss = start + n1 + n2
                                        for p in range( ss + anchor, ss + n3 -anchor+1 ):
                                                if(a2[2],p) in ppL:
                                                        for id in ppL[a2[2],p]:
                                                                #included count at left side
                                                                count[id][nn*3] +=1
                                        for p in range( ss + anchor-1, ss + n3 -anchor):
                                                if(a2[2],p) in ppR:
                                                        for id in ppR[a2[2],p]:
                                                                #included count at left side
                                                                count[id][nn*3+1] +=1

                                start += (n1+n2)

		nn += 1
		fr2.close()
	fr1 = open(countfile)
	fw = open(output,"w")
	for info1 in fr1:
		a1 = info1.strip().split("\t")
		a = intron[a1[0]][0].strip().split()
		fw.write(a1[0] +"\t" + a1[1]+"\t" )
		if(len(intron[a1[0]]) ==1):
			in_p = a1[0].split("_")
			l = int(in_p[2]) -  int(in_p[1])+1
			fw.write(str(l) )
			for i in range(0, num):
				fw.write("\t" + a1[i*6+7]+ "\t"+ a1[i*6+9]+"\t"+a1[i*6+12] )				
		else:
			in1 = [0] * num
			in2 = [0] * num
			in3 = [0] * num
			l = 0
			for j in range(1, len(intron[a1[0]])):
                        	in_p= intron[a1[0]][j].split("_")
                                l += int(in_p[2]) -  int(in_p[1])+1
			fw.write(str(l) )
			if(l ==1):
				for i in range(0, num):
					fw.write("\t0\t0\t0")
				fw.write("\n")
				continue			
			for k in range(0, num):
				for j in range(1, len(intron[a1[0]])):
					if(intron[a1[0]][j] in intron):
						temp = intron[intron[a1[0]][j]][0].strip().split()
						# here is the sum of three column, strand info is not considered here. 
						in1[k]+= int(temp[k*6 +7])
						in2[k]+= int(temp[k*6 +9])
						in3[k]+= int(temp[k*6 +12])
					else:		
						in1[k]+= count[intron[a1[0]][j]][k*3]
						in2[k]+= count[intron[a1[0]][j]][k*3+1]
						in3[k]+= count[intron[a1[0]][j]][k*3+2]
				fw.write("\t"+ str(in1[k])+ "\t"+ str(in2[k])+ "\t"+str(in3[k]))	
		fw.write("\n")
	fw.close()
	exit()

print "this is library with strand info"
intron = dict()
fr = open(countfile)
for info in fr:
	a = info.strip().split("\t")
	intron[a[0]] = [info]
fr.close() 
count = dict()
#the dict to store the intron count,first is the inclusion at left side, second is the inclusion at right side, third is the intron counts
ppL = dict()
#the dict to store left edge  position of intron
ppR = dict()
#the dict to store right  edge  position of intron
pos = dict()
# the dictionary to store the introns in a certain windows
bin = 1000
# set the bin size
samfile = sam.split(",")
num = len(samfile)
fr1 = open(input)
for info1 in fr1:
	a1 = info1.strip().split("\t")
	for key in a1[1:]:
		intron[a1[0]].append(key)
		strand = intron[a1[0]][0].split("\t")[2]
		if(key in count):
			continue
		if(key in intron):
			continue
		count[key]= [0]*num*3
		count[key].append("true")
		count[key].append(strand)
		in_p = key.split("_")
		ppL.setdefault((in_p[0],int(in_p[1])),[]).append(key)
		ppR.setdefault((in_p[0],int(in_p[2])),[]).append(key)
		index_s =  int(in_p[1]) / bin
		index_e =  int(in_p[2]) / bin
		for i in  range(index_s , index_e+1):
			pos.setdefault((in_p[0],i),[]).append(key)


fr1.close()
samfile = sam.split(",")
nn = 0
print "library building complete"
for file in samfile:
	fr2 = open(file)
        for info2 in fr2:
                if (info2[0]== "@"):
                        continue
                strand=''
                a2 = info2.strip().split("\t")
                if(read =="P"):
                        if((int(a2[1])/2 )%2 ==0):
                                continue
                        if(re.search("NH:i:1[0-9]$",info2)):
                                continue
                        if(re.search("NH:i:[^1]$",info2)):
                                continue
                        if(re.search("D|H|S",a2[5])):
                                continue
                        if(re.search("I",a2[5])):
                                continue
                if(read =="S"):
                        if(re.search("NH:i:1[0-9]$",info2)):
                                continue
                        if(re.search("NH:i:[^1]$",info2)):
                                continue
                        if(re.search("D|H|S",a2[5])):
                                continue
                        if(re.search("I",a2[5])):
                                continue
		strand =''
		flag = str(length) +'M' 
                ss = int(a2[1]) /16 % 2
                        
                if ((read =="S") & (lib=="first")):
                        if (ss==0):
                                strand="-"
                        else:
                                strand="+"
                if ((read =="S") & (lib=="second")):
                        if (ss==0):
                                strand="+"
                        else:
                                strand="-"
                
                if ((read =="P") & (lib =="first")):
                        f = int(a2[1]) /64 % 2
                        s = int(a2[1]) /128 % 2
                        if ((ss == 0) & (f ==1)):
                                strand = "-"
                        if ((ss == 0) & (s ==1)):
                                strand = "+"
                        if ((ss == 1) & (f ==1)):
                                strand = "+"
                        if ((ss == 1) & (s ==1)):
                                strand = "-"
                
                if ((read =="P") & (lib =="second")):
                        f = int(a2[1]) /64 % 2
                        s = int(a2[1]) /128 % 2
                        if ((ss == 0) & (f ==1)):
                                strand = "+"
                        if ((ss == 0) & (s ==1)):
                                strand = "-"
                        if ((ss == 1) & (f ==1)):
                                strand = "-"
                        if ((ss == 1) &( s ==1)):
                                strand = "+"
                if (strand==''):
                        print a2
		aa1 = a2[5].split("M")
                aa2 = a2[5].split("N")
                l = len(aa2)
                if(len(aa1)== 2):
                        for p in range(int(a2[3]) + anchor, int(a2[3]) + int(aa1[0]) - anchor +1):
                                if(a2[2],p) in ppL:
                                        for id in ppL[a2[2],p]:
                                                if(re.search(("\%s" %strand),count[id][num*3+1])):
                                                        #inclusion count at left side
                                                        count[id][nn*3] +=1
                                                        continue
                        for p in range(int(a2[3]) + anchor-1, int(a2[3]) + int(aa1[0]) - anchor):
                                if(a2[2],p) in ppR:
                                        for id in ppR[a2[2],p]:
                                                if(re.search(("\%s" %strand),count[id][num*3+1])):
                                                        #inclusion count at right side
                                                        count[id][nn*3+1] +=1 
                                                        continue        
                        index= int(a2[3])/bin
                        if (a2[2], index) in pos:
                                for key in pos[a2[2],index]:
                                        in_p = key.split("_")
                                        if ((int(in_p[1]) <= int(a2[3]) )& (int(in_p[2])>= (int(a2[3]) + length -1))):
                                                if (re.search(("\%s" %strand),count[key][num*3+1])):
                                                        count[key][nn*3+2]+=1
			
       
		start = int(a2[3])
                for i in range(0,l-1):
                        if(re.search("\D",aa2[i].split("M")[0])):
                                continue
                        if(re.search("\D",aa2[i].split("M")[1])):
                                continue
                        if(re.search("\D",aa2[i+1].split("M")[0])):
                                continue
                        n1=int(aa2[i].split("M")[0])
                        n2=int(aa2[i].split("M")[1])
                        n3=int(aa2[i+1].split("M")[0])

		    	if( n1 -anchor+1 >  anchor):
                                ss = start
                                for p in range( ss + anchor, ss + n1 -anchor+1 ):
                                        if(a2[2],p) in ppL:
                                                for id in ppL[a2[2],p]:
                                                        if(re.search(("\%s" %strand),count[id][num*3+1])):
                                                                #included count at left side
                                                                count[id][nn*3] +=1
                                for p in range( ss + anchor-1, ss + n1 -anchor ):
                                        if(a2[2],p) in ppR:
                                                for id in ppR[a2[2],p]:
                                                        if(re.search(("\%s" %strand),count[id][num*3+1])):
                                                                #included count at right side
                                                                count[id][nn*3+1] +=1
                        if( n3 -anchor+1 >  anchor):
                                ss = start + n1 + n2
                                for p in range( ss + anchor, ss + n3 -anchor+1 ):
                                        if(a2[2],p) in ppL:
                                                for id in ppL[a2[2],p]:
                                                        if(re.search(("\%s" %strand),count[id][num*3+1])):
                                                                #included count at left side
                                                                count[id][nn*3] +=1
                                for p in range( ss + anchor-1, ss + n3 -anchor):
                                        if(a2[2],p) in ppR:
                                                for id in ppR[a2[2],p]:
                                                        if(re.search(("\%s" %strand),count[id][num*3+1])):
                                                                #included count at left side
                                                                count[id][nn*3+1] +=1
			
			start += (n1+n2)
	nn += 1
	fr2.close()

fr1 = open(countfile)
fw = open(output,"w")
for info1 in fr1:
	a1 = info1.strip().split("\t")
	a = intron[a1[0]][0].strip().split()
	fw.write(a1[0] +"\t" + a1[1]+"\t" )
	if(len(intron[a1[0]]) ==1):
		in_p = a1[0].split("_")
		l = int(in_p[2]) -  int(in_p[1])+1
		fw.write(str(l) )
		for i in range(0, num):
			fw.write("\t" + a1[i*6+7]+ "\t"+ a1[i*6+9]+"\t"+a1[i*6+12] )                            
	else:
		in1 = [0] * num
		in2 = [0] * num
		in3 = [0] * num
		l = 0
		for j in range(1, len(intron[a1[0]])):
			in_p= intron[a1[0]][j].split("_")
			l += int(in_p[2]) -  int(in_p[1])+1
		fw.write(str(l) )
		if(l ==1):
			for i in range(0, num):
				fw.write("\t0\t0\t0")
			fw.write("\n")
			continue                        
		for k in range(0, num):
			for j in range(1, len(intron[a1[0]])):
				if(intron[a1[0]][j] in intron):
					temp = intron[intron[a1[0]][j]][0].strip().split()
					# here is the sum of three column, strand info is not considered here. 
					in1[k]+= int(temp[k*6 +7])
					in2[k]+= int(temp[k*6 +9])
					in3[k]+= int(temp[k*6 +12])
				else:           
					in1[k]+= count[intron[a1[0]][j]][k*3]
					in2[k]+= count[intron[a1[0]][j]][k*3+1]
					in3[k]+= count[intron[a1[0]][j]][k*3+2]
			fw.write("\t"+ str(in1[k])+ "\t"+ str(in2[k])+ "\t"+str(in3[k]))        
	fw.write("\n")
fw.close()
exit()

