## Representative RNA-seq Scripts

### Build Genome index using STAR

```bash
module load star/2.5.3a
cd ./Genomes/mm9
STAR --runMode genomeGenerate \
--genomeDir ./starIndex \
--genomeFastaFiles ./mm9.fa \
--sjdbOverhang 100 \
--runThreadN 8 \
--sjdbGTFfile ./gencode.mouse.v1.annotation.gtf
```



### Detecting gene expression using STAR

```bash
module load star/2.5.3a

runSTAR(){
## fq1, fq2, outputFile1, outputFile1\SJ.out.tab, outputFile2, Genome4Pass2File, Thread 
STAR \
--genomeDir /rsrch3/scratch/genomic_med/dhao/Genomes/mm9/starIndex \
--readFilesIn $1 $2 \
--outFilterMultimapScoreRange 1 \
--outFilterMultimapNmax 10 \
--outFilterMismatchNoverReadLmax 0.05 \
--outFilterMismatchNoverLmax 0.05 \
--runThreadN $7 \
--alignIntronMax 500000 \
--alignMatesGapMax 1000000 \
--sjdbScore 2 \
--alignSJDBoverhangMin 1 \
--genomeLoad NoSharedMemory \
--readFilesCommand zcat \
--outFilterMatchNminOverLread 0.33 \
--outFilterScoreMinOverLread 0.33 \
--sjdbOverhang 100 \
--outSAMstrandField intronMotif \
--outSAMtype None \
--outSAMmode None \
--outFileNamePrefix $3

STAR \
--runMode genomeGenerate \
--genomeDir $6 \
--genomeFastaFiles /rsrch3/scratch/genomic_med/dhao/Genomes/mm9/mm9.fa \
--sjdbOverhang 100 \
--runThreadN $7 \
--sjdbFileChrStartEnd $4

STAR \
--genomeDir $6 \
--readFilesIn $1 $2 \
--outFilterMultimapScoreRange 1 \
--outFilterMultimapNmax 10 \
--outFilterMismatchNoverReadLmax 0.05 \
--outFilterMismatchNoverLmax 0.05 \
--alignIntronMax 500000 \
--alignMatesGapMax 1000000 \
--sjdbScore 2 \
--alignSJDBoverhangMin 1 \
--genomeLoad NoSharedMemory \
--limitBAMsortRAM 0 \
--readFilesCommand zcat \
--outFilterMatchNminOverLread 0.33 \
--outFilterScoreMinOverLread 0.33 \
--sjdbOverhang 100 \
--runThreadN $7 \
--outSAMstrandField intronMotif \
--outSAMattributes NH HI NM MD AS XS \
--outSAMunmapped Within \
--outSAMtype BAM SortedByCoordinate \
--outSAMheaderHD @HD VN:1.4 \
--outFileNamePrefix $5 \
--outSAMattrRGline ID:id LB:lib SM:sample PL:platform
}

runSTAR $fq1 $fq2 $outputFile1 $outputFile1\$fastqfolder\SJ.out.tab $outputFile2\$fastqfolder $Genome4Pass2File 20
samtools index $bamfile
```

### Counting reads

```bash
htseq-count -m intersection-nonempty -i gene_id -s reverse -r pos -f bam \
--max-reads-in-buffer 50000000 \
${samplename} ./gencode.mouse.v1.annotation.gtf > $HTseq_output/${ouput%Aligned.sortedByCoord.out.bam}.txt
```

### Gene expression analysis

```R
library(edgeR)
geneAnno = read.table("mm9_ProteinCodingGenes.bed", header = FALSE, sep = "\t",stringsAsFactors = F)
names(geneAnno) = c("chr","left","right","gid","name","strand")# only protein coding genes
geneAnno <- geneAnno[geneAnno$chr!="chrM",] #(MT genes deleted)
rownames(geneAnno) <- geneAnno$gid
############################################################################################
# Load gene expression data
header = read.table("./HTSeq_counts/sample1.txt", header = FALSE,stringsAsFactors = F)
header <- header$V1
iOrd <- list.files("./HTSeq_counts",full.names = T)
ReadCounts <- NULL
for (i in iOrd) {
  x <- read.table(i, header = FALSE,stringsAsFactors = F)
  ReadCounts <- cbind(ReadCounts,x[,2])
}
rownames(ReadCounts) <- header
ReadCounts <- ReadCounts[1:(nrow(ReadCounts)-5),]
iOrd <- basename(iOrd);iOrd<- gsub(".txt","",iOrd)
colnames(ReadCounts) <- paste0("s",iOrd)
ReadCounts <- ReadCounts[geneAnno$gid,] # only keep protein coding genes
save(ReadCounts,file="./ReadCounts.rda")

EdgeR <- function(ReadCounts, group1, group2) {
  # Create design matrix
  condition <- c(rep("group1",length(group1)),rep("group2",length(group2)))
  condition <- factor(condition,levels = c("group1","group2"))
  design = model.matrix(~condition)
  # load gene expression data:
  ReadCounts <- ReadCounts[,c(group1,group2)]
  geneList = DGEList(counts=ReadCounts, genes=rownames(ReadCounts))
  # filter lowly expressed genes/transcripts and recompute the library sizes:
  #geneList = geneList[rowSums(ReadCounts)>10, , keep.lib.sizes=FALSE] # at least 10 reads on all samples
  # Apply TMM normalization:
  geneList = calcNormFactors(geneList, method="TMM") #"TMMwzp"performs better for data with a high proportion of zeros
  # Estimating dispersion
  geneList = estimateDisp(geneList, design,robust=TRUE)
  # Differential expression:
  fit <- glmFit(geneList, design)
  lrt <- glmLRT(fit)
  # Gnerate tabular output
  geneDE = topTags(lrt, n = nrow(lrt),adjust.method ="fdr")$table
  geneDE$LR <- ifelse(geneDE$logFC>0,geneDE$LR,-1*geneDE$LR)
  return(geneDE)
}

group2 <- c("s3a1KO_CTX_1","s3a1KO_CTX_2","s3a1KO_CTX_3")
group1 <- c("sWT_CTX_1","sWT_CTX_2","sWT_CTX_3")
KO3a1.DEG <- EdgeR(ReadCounts,group1,group2)
KO3a1.DEG$Symbol <- geneAnno[rownames(KO3a1.DEG),"name"]
KO3a1.DEGs <- KO3a1.DEG[KO3a1.DEG$PValue < p_cutoff,]

# location of DEGs along the BED file sorted by Dnmt3a binding:
BED <- fread("./byDnmt3a.bed")
BED <- data.frame(BED,stringsAsFactors = F)
rownames(BED) <- BED$name
KO3a1.DEGs <- KO3a1.DEGs[is.element(KO3a1.DEGs$genes,BED$name),]
BED$KO3a1DEG <- NA;BED[KO3a1.DEGs$genes,"KO3a1DEG"] <- 1
BED$type <- NA;BED[KO3a1.DEGs$genes[KO3a1.DEGs$logFC>0],"type"] <- "upDEG"
BED[KO3a1.DEGs$genes[KO3a1.DEGs$logFC<0],"type"] <- "dnDEG"
ggplot() + 
geom_segment(data=BED,mapping=aes(x=seq_along(BED$type), y=0, xend=seq_along(BED$type),yend=KO3a1DEG,color=type),size=0.1)+
theme(axis.title.x=element_blank(),axis.text.x=element_blank(),axis.ticks.x=element_blank(),
        panel.border = element_blank(), panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black")) + 
theme_bw()
```

### Correlation between methylation change and expression change

```bash
# obtain methylation value at gene body
multiBigwigSummary BED-file --binSize 100000 \
 -b sample_01.bw sample_02.bw sample_03.bw sample_04.bw \
 --labels WT_1 WT_2 3a1KO_1 3a1KO_2
 -out scores_per_bin.npz \
 --outRawCounts scores_per_GeneBody_CpG.tab -p 5 \
 --BED ProteinCodingGenes.bed
```

```R
library(data.table);library(ggplot2)
WGBS <- fread("scores_per_GeneBody_CpG.tab",header = T)

Sliding <- function(TotalNumber, binSize, SlidingSize) {
  bins <- list()
  i=1;iOrd_start=1;iOrd_end=binSize
  while(iOrd_end < TotalNumber){
    bins[[i]] <- iOrd_start:iOrd_end
    i = i+1
    iOrd_start = iOrd_start+SlidingSize
    iOrd_end <- ifelse((iOrd_end + SlidingSize) < TotalNumber,iOrd_end + SlidingSize,TotalNumber)
  }
  bins[[i]] <- iOrd_start:iOrd_end
  return(bins)
}

KO3a1.WGBS <- WGBS[WGBS$gene %in% KO3a1.DEG$genes,]
KO3a1.WGBS$meanDif <- rowMeans(KO3a1.WGBS[,c("3a1KO_NeuNuclei_Me_1","3a1KO_NeuNuclei_Me_2")]) - rowMeans(KO3a1.WGBS[,c("WT_NeuNuclei_Me_1","WT_NeuNuclei_Me_2")])
  KO3a1.WGBS$meanDif <- KO3a1.WGBS$meanDif*100 # % percentages
  KO3a1.WGBS <- KO3a1.WGBS[order(KO3a1.WGBS$meanDif,decreasing = F),]
  KO3a1.WGBS <- KO3a1.WGBS[!is.na(KO3a1.WGBS$meanDif),]
  KO3a1.WGBS$logFC <- KO3a1.DEG[KO3a1.WGBS$gene,"logFC"]
  bins <- Sliding(nrow(KO3a1.WGBS),100,20)
  KO3a1.cor <- data.frame(iOrd=1:length(bins),logFC=NA,percent_methy=NA)
  for (i in 1:length(bins)) {
    KO3a1.cor$iOrd[i] <- i
    KO3a1.cor$logFC[i] <- mean(KO3a1.WGBS$logFC[bins[[i]]])
    KO3a1.cor$percent_methy[i] <- mean(KO3a1.WGBS$meanDif[bins[[i]]])
  }
  
  p1 <- ggplot(KO3a1.cor, aes(x = percent_methy, y = logFC)) + geom_point(color = "#00AFBB",size=2) + 
    geom_smooth(method = lm, se = FALSE) + theme_classic()+
    ylab("average LogFC (3a1Ko vs WT)") + xlab("%mCpG (3a1KO-WT)")+
```



## Softwares

STAR (v2.5.3b) 

HTSeq-count (v0.11.0) 

R package EdgeR (v3.26.8), ggplot2 (v3.3.5) 