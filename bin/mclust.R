args=commandArgs(TRUE)
print(args)
if(length(args)) {
	FileName=args[1]
}else{ 
	FileName="infile"
}
library(mclust)
data<-read.table(FileName,head=T,sep="\t",check.names=F)
dd <- dim(data)
#head(data)
samplenum <- dd[2]-10
addcol <- dd[2]-1
matrix_output <- matrix(nrow = dd[1], ncol = addcol )
colnumber <- dd[2]-2
for (i in 1:dd[1]){
	if (i %% 1000 == 0) print(i)
	mclustinput <- as.numeric(data[i,9:colnumber])
	if(length(table(mclustinput)) == 1){
		classfication <- rep(1,samplenum)
		means <- rep(mclustinput[1],3)
		variances <- rep(0,3)
	}else{
		mod<-densityMclust(mclustinput,G=1:3)###set components number range
		classfication <- as.vector(summary(mod,classification = T)$classification)
		mean <- as.vector(summary(mod,parameters = T)$mean)
		variance <- as.vector(summary(mod,parameters = T)$variance)
		means <- c(mean[1], mean[2], mean[3])
		variances <- c(variance[1], variance[2], variance[3])
	}
	matrix_output[i,] <- c(classfication, means, variances, sum(classfication == 1), sum(classfication == 2), sum(classfication == 3))
}
output <- cbind(data,matrix_output)
write.table(output,file="Genotype_Unsupervised",col.names =T, row.names = F,quote = F,sep = "\t")
