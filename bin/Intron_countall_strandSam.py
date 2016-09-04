#!/bin/python
import getopt,copy,re,os,sys,logging,time,datetime;
options, args = getopt.getopt(sys.argv[1:], 'o:',['gtf=','anchor=','lib=','read=','length=','sam=','output=','Total='])
gtf='';
read='';
sam='';
lib='';
length =0
anchor =0
output='';
Total='';
for opt, arg in options:
	if opt in ('-o','--output'):
		output = arg
	elif opt in ('--gtf'):
		gtf = arg
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
	elif opt in ('--Total'):
                Total= arg
if (not gtf or not read or not sam or not output or not length or not anchor or not Total):
	print "Not enough parameters!"
        print "Program : ", sys.argv[0]
        print "          A python program to count the reads for retained intron events for varities of junction from a series of sam file."
        print "Usage :", sys.argv[0], " --gtf: the intron gtf file;"
        print "Usage :", sys.argv[0], " --length:the length of reads;"
        print "Usage :", sys.argv[0], " --anchor:the anchor length of the read;"
        print "Usage :", sys.argv[0], " --sam: the sam file,multiple sam file seperated by commas;"
        print "Usage :", sys.argv[0], " --lib: the library type;"
        print "Usage :", sys.argv[0], " --read: The sequencing strategy of producing reads with choices of P/S;"
        print "Usage :", sys.argv[0], ' --output: intron_id, gene_id,strand,chr,start,end,5SS inclusion counts,5SS skipping counts,3SS includion counts,3SS skipping counts,skipping counts,intron counts;'
	print "Usage :", sys.argv[0], " --Total: the file store the total uniquely mapped reads."
	print datetime.datetime.now()
	print "Author  : Shaofang Li"
	print "Contact : sfli001@gmail.com"
	sys.exit()
if(lib=="unstrand"):
	print "This is the library without strand information"

	fr1 = open(gtf)
	count = dict()
	#the dict to store the intron count,first is the inclusion at left side, second is the skipped count at smaller side, third inclusion at right side, fourth skipping at right side, fifith skipping counts, sixth intron counts, at the end of the three column store gene_id, gene_strand, and clean introns based on bam files
	ppL = dict()
	#the dict to store left edge  position of intron
	ppR = dict()
	#the dict to store right  edge  position of intron
	pp = dict()
	#the dict to store the two edge position of intron
	pos = dict()
	# the dictionary to store the introns in a certain windows
	bin =1000
	# set the bin size
	samfile = sam.split(",")
	num = len(samfile)
	for info1 in fr1:
		a1 = info1.strip().split("\t")
		key = "%s_%s_%s" % (a1[0],a1[3],a1[4])
		gene_id = re.sub('.*gene_id "|\".*','',a1[8])
		if(count.has_key(key)):
			gg = count[key][6*num].split(",")
			gg_l = "true"
			for g_id in gg:
				if(gene_id == g_id):
					gg_l = "falase"
			if(gg_l =="true"):
				count[key][6*num]+=","
				count[key][6*num]+=gene_id
				count[key][6*num+1]+=","
				count[key][6*num+1]+=a1[6]

		else:
			count[key]= [0]*num*6
			count[key].append(gene_id)
			count[key].append(a1[6])
			count[key].append("true")
			count[key].append("true")
			ppL.setdefault((a1[0],int(a1[3])),[]).append(key)
                	ppR.setdefault((a1[0],int(a1[4])),[]).append(key)
			pp[a1[0],int(a1[3]),int(a1[4])] = key
			index_s =  int(a1[3]) / bin
               		index_e =  int(a1[4]) / bin
                	for i in  range(index_s , index_e+1):
                        	pos.setdefault((a1[0],i),[]).append(key)
	fr1.close()
	samfile = sam.split(",")
	nn = 0
	Lib_T = [0] * len(samfile)
	for file in samfile:
		fr2 = open(file)
		T = 0
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
                     	T+=1
			aa1 = a2[5].split("M")
                        aa2 = a2[5].split("N")
                        l = len(aa2)
			if(len(aa1)== 2):
				for p in range(int(a2[3]) + anchor, int(a2[3]) + int(aa1[0]) - anchor +1):
					if(a2[2],p) in ppL:
						for id in ppL[a2[2],p]:
							#inclusion count at left side
							count[id][nn*6] +=1
							continue
				for p in range(int(a2[3]) + anchor-1, int(a2[3]) + int(aa1[0]) - anchor):
					if(a2[2],p) in ppR:
                                                for id in ppR[a2[2],p]:
                                                        #inclusion count at right side
                                                        count[id][nn*6+2] +=1 
                                                        continue	
				index= int(a2[3])/bin
                                if (a2[2], index) in pos:
                                        for key in pos[a2[2],index]:
						in_p = key.split("_")
                                                if ((int(in_p[1]) <= int(a2[3]) )& (int(in_p[2])>= (int(a2[3]) + length -1))):
                                                        count[key][nn*6+5]+=1
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
				if (n1 >= anchor and n3 >= anchor):
					index1 = (start+ n1-1)/ bin
                                        index2 = (start+ n1+n2)/bin
					if (a2[2], index1) in pos:
                                                for key in pos[a2[2],index1]:
                                                        in_p = key.split("_")
                                                        if ((int(in_p[1]) < start +n1 )& (int(in_p[2])> start +n1)):
                                                                count[key][num*6+2] ="false"
                                        if( a2[2], index2) in pos:
                                                for key in pos[ a2[2],index2]:
                                                        in_p = key.split("_")   
                                                        if ((int(in_p[1]) < start +n1+n2 -1)& (int(in_p[2])> start +n1+n2-1)):
                                                                count[key][num*6+2] ="false"		

					if (a2[2], start+n1) in ppL:
						for id in ppL[a2[2],start +n1]:
							#skipped count at left side
							count[id][nn*6+1] +=1
					if (a2[2], start+n1+n2-1) in ppR:
                                                for id in ppR[a2[2],start +n1+n2 -1]:
                                                        #skipped count at right side
                                                        count[id][nn*6+3] +=1
					if (a2[2],start+n1, start+n1+n2-1) in pp:
						#skipping count at both side
        					count[pp[a2[2],start+n1,start +n1+n2 -1]][nn*6+4] +=1	

				if( n1 -anchor+1 >  anchor):
					ss = start
					for p in range( ss + anchor, ss + n1 -anchor+1 ):
						if(a2[2],p) in ppL:
							for id in ppL[a2[2],p]:
								#included count at left side
								count[id][nn*6] +=1
					for p in range( ss + anchor-1, ss + n1 -anchor ):
                                                if(a2[2],p) in ppR:
                                                        for id in ppR[a2[2],p]:
								#included count at right side
                                                                count[id][nn*6+2] +=1
				if( n3 -anchor+1 >  anchor):
					ss = start + n1 + n2
					for p in range( ss + anchor, ss + n3 -anchor+1 ):
						if(a2[2],p) in ppL:
							for id in ppL[a2[2],p]:
								#included count at left side
								count[id][nn*6] +=1
					for p in range( ss + anchor-1, ss + n3 -anchor):
						if(a2[2],p) in ppR:
                                                        for id in ppR[a2[2],p]:
                                                                #included count at left side
                                                                count[id][nn*6+2] +=1

				start += (n1+n2)
		Lib_T[nn] = T
		nn += 1
		fr2.close()
	fr1 = open(gtf)
	fw = open(output,"w")
	for info1 in fr1:
		a1 = info1.strip().split("\t")
		key = "%s_%s_%s" % (a1[0],a1[3],a1[4])
		if(count[key][num*6+3]=="true"):
			if(re.search('\+',count[key][num*6+1])):
				fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (key, count[key][num*6],count[key][num*6+1],count[key][num*6+2],a1[0],a1[3],a1[4],"\t".join(str(x) for x in count[key][0:num*6])))
			else:
				fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (key, count[key][num*6],count[key][num*6+1],count[key][num*6+2],a1[0],a1[3],a1[4]))
				for i in range(0, num):
					fw.write("\t%s\t%s\t%s\t%s\t%s\t%s"%(count[key][i*6+2],count[key][i*6+3],count[key][i*6],count[key][i*6+1],count[key][i*6+4],count[key][i*6+5]))
				fw.write("\n") 
			count[key][num*6+3]="false"
	fw.close()
	fw = open(Total,"w")
        fw.write("\t".join(str(x) for x in Lib_T))
        fw.close()
	exit()

print "this is library with strand info"
fr1 = open(gtf)
count = dict()
#the dict to store the intron count,first is the inclusion at left side, second is the skipped count at smalller side, third inclusion at right side, fourth skipping at right side
ppL = dict()
#the dict to store left edge  position of intron
ppR = dict()
#the dict to store right  edge  position of intron
pp = dict()
#the dict to store the two edge position of intron
pos = dict()
# the dictionary to store the introns in a certain windows
bin =1000
# set the bin size

samfile = sam.split(",")
num = len(samfile)
for info1 in fr1:
	a1 = info1.strip().split("\t")
	key = "%s_%s_%s" % (a1[0],a1[3],a1[4])
        gene_id = re.sub('.*gene_id "|\".*','',a1[8])
	if(count.has_key(key)):
		gg = count[key][6*num].split(",")
		gg_l = "true"
		for g_id in gg:
			if(gene_id == g_id):
				gg_l = "falase"
		if(gg_l =="true"):
			count[key][6*num]+=","
			count[key][6*num]+=gene_id
			count[key][6*num+1]+=","
			count[key][6*num+1]+=a1[6]

	else:
		count[key]= [0]*num*6
		count[key].append(gene_id)
		count[key].append(a1[6])
		count[key].append("true")
		count[key].append("true")
		ppL.setdefault((a1[0],int(a1[3])),[]).append(key)
		ppR.setdefault((a1[0],int(a1[4])),[]).append(key)
		pp[a1[0],int(a1[3]),int(a1[4])] = key
		index_s =  int(a1[3]) / bin
		index_e =  int(a1[4]) / bin
                for i in  range(index_s , index_e+1):
                	pos.setdefault((a1[0],i),[]).append(key)

fr1.close()
samfile = sam.split(",")
nn = 0
Lib_T = [0] * len(samfile)
for file in samfile:
        fr2 = open(file)
	T = 0
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

             	T +=1
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
						if(re.search(("\%s" %strand),count[id][num*6+1])):
							#inclusion count at left side
							count[id][nn*6] +=1
							continue
			for p in range(int(a2[3]) + anchor-1, int(a2[3]) + int(aa1[0]) - anchor):
				if(a2[2],p) in ppR:
					for id in ppR[a2[2],p]:
						if(re.search(("\%s" %strand),count[id][num*6+1])):
							#inclusion count at right side
							count[id][nn*6+2] +=1 
							continue        
			index= int(a2[3])/bin
			if (a2[2], index) in pos:
				for key in pos[a2[2],index]:
					in_p = key.split("_")
					if ((int(in_p[1]) <= int(a2[3]) )& (int(in_p[2])>= (int(a2[3]) + length -1))):
						if (re.search(("\%s" %strand),count[key][num*6+1])):
							count[key][nn*6+5]+=1
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
			if (n1 >= anchor and n3 >= anchor):
				index1 = (start+ n1-1)/ bin
				index2 = (start+ n1+n2)/bin
				if (a2[2], index1) in pos:
					for key in pos[a2[2],index1]:
						in_p = key.split("_")
						if ((int(in_p[1]) < start +n1)& (int(in_p[2])>start +n1)):
							if (re.search(("\%s" %strand),count[key][num*6+1])):
								count[key][num*6+2] ="false"
				if( a2[2], index2) in pos:
					for key in pos[ a2[2],index2]:
						in_p = key.split("_")
						if ((int(in_p[1]) < start +n1+n2 -1)& (int(in_p[2])> start +n1+n2-1)):
							if (re.search(("\%s" %strand),count[key][num*6+1])):
								count[key][num*6+2] ="false"
				if (a2[2], start+n1) in ppL:
					for id in ppL[a2[2],start +n1]:
						if(re.search(("\%s" %strand),count[id][num*6+1])):
							#skipped count at left side
							count[id][nn*6+1] +=1
				if (a2[2], start+n1+n2-1) in ppR:
					for id in ppR[a2[2],start +n1+n2 -1]:
						if(re.search(("\%s" %strand),count[id][num*6+1])):
							#skipped count at right side
							count[id][nn*6+3] +=1
				if (a2[2],start+n1, start+n1+n2-1) in pp:
					if(re.search(("\%s" %strand),count[pp[a2[2],start+n1,start +n1+n2 -1]][num*6+1])):
						#skipping count at both side
						count[pp[a2[2],start+n1,start +n1+n2 -1]][nn*6+4] +=1   
		    	if( n1 -anchor+1 >  anchor):
				ss = start
				for p in range( ss + anchor, ss + n1 -anchor+1 ):
					if(a2[2],p) in ppL:
						for id in ppL[a2[2],p]:
							if(re.search(("\%s" %strand),count[id][num*6+1])):
								#included count at left side
								count[id][nn*6] +=1
				for p in range( ss + anchor-1, ss + n1 -anchor ):
					if(a2[2],p) in ppR:
						for id in ppR[a2[2],p]:
							if(re.search(("\%s" %strand),count[id][num*6+1])):
								#included count at right side
								count[id][nn*6+2] +=1
			if( n3 -anchor+1 >  anchor):
				ss = start + n1 + n2
				for p in range( ss + anchor, ss + n3 -anchor+1 ):
					if(a2[2],p) in ppL:
						for id in ppL[a2[2],p]:
							if(re.search(("\%s" %strand),count[id][num*6+1])):
								#included count at left side
								count[id][nn*6] +=1
				for p in range( ss + anchor-1, ss + n3 -anchor):
					if(a2[2],p) in ppR:
						for id in ppR[a2[2],p]:
							if(re.search(("\%s" %strand),count[id][num*6+1])):
								#included count at left side
								count[id][nn*6+2] +=1

			start += (n1+n2)
	Lib_T[nn] = T
	nn += 1
	fr2.close()

fr1 = open(gtf)
fw = open(output,"w")
for info1 in fr1:
	a1 = info1.strip().split("\t")
	key = "%s_%s_%s" % (a1[0],a1[3],a1[4])
	if(count[key][num*6+3]=="true"):
		if(re.search('\+',count[key][num*6+1])):
			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" % (key, count[key][num*6],count[key][num*6+1],count[key][num*6+2],a1[0],a1[3],a1[4],"\t".join(str(x) for x in count[key][0:num*6])))
		else:
			fw.write("%s\t%s\t%s\t%s\t%s\t%s\t%s" % (key, count[key][num*6],count[key][num*6+1],count[key][num*6+2],a1[0],a1[3],a1[4]))
			for i in range(0, num):
				fw.write("\t%s\t%s\t%s\t%s\t%s\t%s"%(count[key][i*6+2],count[key][i*6+3],count[key][i*6],count[key][i*6+1],count[key][i*6+4],count[key][i*6+5]))
			fw.write("\n") 
		count[key][num*6+3]="false"
fw.close()
fw = open(Total,"w")
fw.write("\t".join(str(x) for x in Lib_T))
fw.close()

