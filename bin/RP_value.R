args <- commandArgs(trailingOnly = TRUE)
#args[1] use this as the input Diff_compare1_intron_PI.txt
#args[2], the time of permutation perform
#args[3], the output files
#args[4], the output files
#args = c("/u/nobackup/yxing/NOBACKUP/shaofang/shRNA/K562/new/squid_0.1/SQUID_U2AF2_K562_ENCSR904CJQ/Result/Diff_IR_1_intron_PI.txt","100")
times = as.numeric(args[2])
a = read.table(args[1],header = TRUE)
data = a[complete.cases(a),]
data = data[,c(1,15,24)]


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
write.table(data,args[3],row.names = FALSE,quote = FALSE, sep="\t")
write.table(re,args[4],row.names = FALSE,quote = FALSE, sep="\t")
