#BiocManager::install("Rsubread")  #Install package of Rsubread
library("Rsubread")   #load Rsubread
library(edgeR)  # load edgeR
#devtools::install_github("zhangyuqing/sva-devel")
library(sva) 
setwd("d:/37662/deg/")
samples=read.csv("../samples37662.csv",header=T)
samples

ref <- "d:/37662/final_ref/wt_37662.fa"   # 参考基因组路径及文件
#建立索引，两个参数，basename即索引名字，自己取名，会在后续比对时用到。reference是文件路径
buildindex(basename="37662wt",reference=ref) 
#序列的路径
bgireads_file="D:\\37662\\bgi_transcriptome_37662\\01.cleanData/"
aptreads_file="D:\\37662\\multi-omics\\P20191202270_transcriptome_37662/"
samples$file1[1:3]=paste0(bgireads_file,samples$file1[1:3])
samples$file2[1:3]=paste0(bgireads_file,samples$file2[1:3])
samples$file1[4:9]=paste0(aptreads_file,samples$file1[4:9])
samples$file2[4:9]=paste0(aptreads_file,samples$file2[4:9])
samples

# 比对开始
for(i in 1:9){
 align(index="37662wt",readfile1=samples$file1[i],readfile2=samples$file2[i],output_file=paste0(samples$samples[i],".bam"),nthreads = 10)
}
fc_se=list()
for(i in 1:9){
  fc_se[[i]] <- featureCounts(paste0(samples$samples[i],".bam"),annot.ext="../pgap/wt37662.gff3",isGTFAnnotationFile=T,
                            GTF.featureType = "gene",GTF.attrType = "Name",isPairedEnd=TRUE,nthreads = 10)
}

groups= samples$group
groups=factor(groups,levels=c("wt","rm1","rm2"),order=T)
group <- factor(c(1,2,3,1,2,3,1,2,3))
counts=fc_se[[1]]$counts
for(i in 2:9){
counts=cbind(counts,fc_se[[i]]$counts)
}
head(counts)


# To see if there are batch effects
y <- DGEList(counts=counts, group=group)
keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
plotMDS(y, col=rep(1:3, each=3))

### it did exist batch effect, remove it by ComBat-seq
#  devtools::install_github("zhangyuqing/sva-devel")  
library("sva")
adjusted <- ComBat_seq(counts, batch=samples$batch, group=group)
y <- DGEList(counts=adjusted, group=group)
keep <- filterByExpr(y)
table(keep)
y <- y[keep, , keep.lib.sizes=FALSE]
y <- calcNormFactors(y)
y$samples
### plot again to see if batch effect removed
plotMDS(y, col=rep(1:3, each=3))

design <- model.matrix(~group)
y <- estimateDisp(y, design)
fit <- glmQLFit(y, design)
qlf <- glmQLFTest(fit, coef=2:3)
ab=topTags(qlf,n=100)
write.csv(ab,file="DEgenes.csv")
cpm=counts[rownames(ab),]
write.csv(cpm,file="DEgenes_cpm.csv")

