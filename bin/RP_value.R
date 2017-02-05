args <- commandArgs(trailingOnly = TRUE)
#args[1] use this as the input Diff_compare1_intron_PI.txt
#args[2], the time of permutation perform
#args[3], the minimum intron length,default is 100
#args[4], the file with the FPKM value of  intron overlapping transcripts
#args[5], the minimum mean FPKM value of intron overlapping transcripts in sample1 and sample2,default is 5
#args[6], the order of sample1 in the comparison file
#args[7], the order of sample2 in the comparison file
#args[8], the minimum mean counts of  inclusion and skipping counts in sample1 and sample2, default is 20
#args[9], the output files of rank product test
#args[10], the final output file in the result folder
#args[11], the cutoff of delta to output differential spliced introns.The default is 0.05
#args[12], the cutoff of combined FDR to output differential spliced introns.The default is 0.05
#args[13], the final output file in the result folder showing the introns with increased PI in sample2
#args[14], the final output file in the result folder showing the introns with decreased PI in sample2
times = as.numeric(args[2])

a = read.table(args[1],header = TRUE)
a= a[complete.cases(a),]
##filter based on intron length
a = a[a[,6]-a[,5] + 1 >= as.numeric(args[3]),]

##filter based on FPKM value
G = read.table(args[4],header = TRUE)
ss1 =  as.numeric(unlist(strsplit(args[6],",")))+1
ss2 =  as.numeric(unlist(strsplit(args[7],",")))+1
G =  G[rowSums(G[,ss1]) > length(ss1)*as.numeric(args[8]) & rowSums(G[,ss2]) > length(ss2)*as.numeric(args[8]),]

a = a[a[,1] %in%G[,1],]
a= a[complete.cases(a),]
##filter based on sum of inclusion and excluion counts

I_s1 = unlist(lapply(a[,9], function(x) mean(as.numeric(unlist(strsplit(as.character(x),","))))))
S_s1 = unlist(lapply(a[,10], function(x) mean(as.numeric(unlist(strsplit(as.character(x),","))))))

I_s2 = unlist(lapply(a[,11], function(x) mean(as.numeric(unlist(strsplit(as.character(x),","))))))
S_s2 = unlist(lapply(a[,12], function(x) mean(as.numeric(unlist(strsplit(as.character(x),","))))))
a = data.frame(a,I_s1,S_s1,I_s2,S_s2 )
a = a[a[,29] + a[,30] > args[8] & a[,31] + a[,32]> args[8],]
a = a[,1:28]
data = a[,c(1,15,24)]


data = data[order(data[,2]),]
data = data.frame(data, 1:length(data[,1]))
data = data[order(data[,3]),]
data = data.frame(data, 1:length(data[,1]))
data =  data.frame(data,RP= as.integer(sqrt(data[,4])*sqrt(data[,5])))
data = data[order(data[,6]),]
data = data.frame(data, 1:length(data[,1]))
len = length(data[,1])
rank = sqrt(1:len)
per = numeric(0)
for(i in 1:times)
{
	temp = data.frame(sample(rank,len),sample(rank,len))
	rp = apply(temp,1, prod)
	per = c(rp, per)
} 
per = as.integer(per[order(per)])
c = numeric(len)
index1 = 1
index2 = times
for (i in 1:len)
{

	temp = per[index1:index2]
	while(length(temp[temp< data[i,6]])==times)
	{
		index1 = index1 + times
		index2 = index2 + times
		temp = per[index1:index2]
	} 	 
	c[i] = index1-1 + length(temp[temp< data[i,6]])
}

data = data.frame(data,c)

pfp = numeric(len)

for (j in 1:len)
{
pfp[j] = data[j,8]/times/data[j,7]
}
pfp[pfp>1]=1
data = data.frame(data,pfp)
data = merge(data,a[,c(1,19,28)],by.x = 1, by.y = 1, all.x = TRUE, all.y = FALSE)
re = merge(a, data[,c(1,9)],by.x =1 ,by.y = 1, all = TRUE)
colnames(re) = c(colnames(a),"Combined_FDR")
colnames(data) = c("Intron_id", "PValue_rMATS", "PValue_DEXSeq","rank_rMATS", "rank_DEXSeq","RP","rank_RP","c","pfp","Diff_PI_Junction","Diff_PI_Density")
write.table(data,args[9],row.names = FALSE,quote = FALSE, sep="\t")
write.table(re,args[10],row.names = FALSE,quote = FALSE, sep="\t")
re1 = re[re[,19] < -as.numeric(args[11]) & re[,28] < -as.numeric(args[11]) & re[,29] < as.numeric(args[12]),]
re2 = re[re[,19] > as.numeric(args[11]) & re[,28] > as.numeric(args[11]) & re[,29] < as.numeric(args[12]),]

write.table(re1,args[13],row.names = FALSE,quote = FALSE, sep="\t")
write.table(re2,args[14],row.names = FALSE,quote = FALSE, sep="\t")

