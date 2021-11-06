## Representative WGBS Scripts

### mm9 genome preparation

```bash
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm9/bigZips/chromFa.tar.gz
tar -xvzf chromFa.tar.gz
rm -rf *_random.fa
cat *.fa > mm9.fa
# index the concatenated .fa file using `samtools` and bwa
module load samtools
samtools faidx mm9.fa
module load bwa
bwa index mm9.fa
```

### Mapping WGBS data to mm9 using bsmap

```bash
module load gcc/7.2.0
module load samtools/0.1.19
mkdir bsmap_out
./software/bsmap-2.90/bsmap -a sample1_R1.fastq -b sample1_R2.fastq -d \
./Genomes/mm9/mm9.fa \
-p 6 -o ./bsmap_out/sample1.bam
samtools sort -@ 5 ./bsmap_out/sample1bam ./bsmap_out/sample1_sorted
```

### Estimate CT conversion rate

```bash
module load moabs
mcall -m $samplename --statsOnly 1 -r ./Genomes/mm9/mm9.fa -p 8 --sampleName ${samplename%_sorted.bam}_bisulfite_conversion
```

### Estimate methylation ratio of CpG sites

```bash
module load samtools/0.1.19

python ./software/bsmap-2.90/methratio.py -d ./Genomes/mm9/mm9.fa -x CG -i no-action -u -g -m 5 -p -r -w ./MethyRatio/${samplename%_sorted.bam}.wig -o ./MethyRatio/${samplename%_sorted.bam}.methyRatio $bamPath/$samplename
```

### De Novo identification of DMR

```bash
#1)formatting data for metilene:
for i in `ls *methyRatio`
do
tail $i -n+2 | awk -F '\t' '{print $1"\t"$2"\t"$2+1"\t"$5}' > $i.txt
done

#2) combine meth Ratio of samples to compare:
module load bedtools
bedtools unionbedg -i sample1.methyRatio.txt sample2.methyRatio.txt -filler - -header > group1.txt
awk -F '\t' '{print $1"\t"$2"\t"$4"\t"$5"\t"$6"\t"$7}' group1.txt > group_1.txt # remove 'end'
sed -i '1 s/^.*/chr\tpos\tc1_1\tc1_2\tt1_1\tt1_2/' group_1.txt # add headers

#3) De novo DMR by metilene
./software/metilene_v0.2-8/metilene_linux64 \
-t 8 -f 1 -c 2 -d 0.1 -a c1 -b t1 -X 1 -Y 1 group_1.txt > ./DMR/group_1_DMR.txt
# add headers to the file
sed -i '1s/^/chr\tstart\tstop\tq-value\tmean difference: mean g1 - mean g2\t#CpGs\tp (MWU)\tp (2D KS)\tmean g1\tmean g2\n/' ./DMR/group_1_DMR.txt
```

### Annotation of DMR using HOMER

```bash
module load homer
annotatePeaks.pl group_1_DMR.bed mm9 -annStats group_1_DMR_status.txt > group_1_DMR_annotated.bed
```

### Identification of UMR (Under-methylated regions)

```bash
module load python/3.7.3-anaconda
./iBSTools_v1.3.0/towig --head --depth 5 --ratio_total 5,6 -i sample1.methyRatio -o ./wigs -n sample1
./iBSTools_v1.3.0/pattern -i ./wigs/sample1.wig -o ./wigs/ -n sample1
```

### Density distribution of CpG methylation

```bash
module load deeptools

# convert wig files to bigwig
./software/wigToBigWig.dms ./MethyRatio/${samplename%_sorted.bam}.wig ./Genomes/mm9/mm9.len ./MethyRatio/${samplename%_sorted.bam}.bw

computeMatrix scale-regions -S ${samplename}.bw\
-R ./protein-coding.bed -b 2000 -a 2000 --outFileName ./matrix.gz \
--regionBodyLength 5000

#plot profile:
plotProfile -m ./deeptool_result/matrix.gz -out ./density_Genes_perGroup.pdf --plotTitle "" --numPlotsPerRow 3 --perGroup
```

### Overlap between DMR and gene body

```R
library(EnrichedHeatmap)
library(RColorBrewer)
library(GenomicRanges)
library(circlize)

Gene <- fread("protein-coding-genes.bed")
Gene <- data.frame(Gene,stringsAsFactors = F)
Gene <- GRanges(seqnames = Gene[,1],
               ranges = IRanges(start = Gene$start, end = Gene$end, names = Gene$name),
               strand = Gene$strand,
               deepTools_group = Gene$deepTools_group,
               orders=1:nrow(Gene))

# heatmap of DMR around gene TSS_ordered by Dnmt3a signal:
DMR <- fread("group_1_DMR.txt")
DMR <- data.frame(DMR,stringsAsFactors = F)
#DMR_3a1KO <- DMR_3a1KO[DMR_3a1KO$q.value<0.1,]
DMR <- DMR_3a1KO[DMR$p..MWU.<0.01,]

DMR <- GRanges(seqnames = DMR$chr,
                     ranges = IRanges(start = DMR$start, end = DMR$stop),
                     strand = "*",
                     dif = DMR$mean.difference..mean.g1...mean.g2)
mat_DMR = normalizeToMatrix(DMR, Gene, mean_mode = "absolute",extend = 5000,target_ratio=0.5)
pdf("./heatmap of DMR across Genes.pdf",width = 3,height = 6)
print(EnrichedHeatmap(mat_DMR, col = c("white", "red"), name = "DMR",row_order=NULL,
                      top_annotation = HeatmapAnnotation(enriched = anno_enriched(ylim=c(0,0.22))),use_raster=F))
dev.off()
```

### Methylation changes across gene body

```bash
module load deeptools/3.1.3
computeMatrix scale-regions -S sample_01.bw sample_02.bw sample_03.bw sample_04.bw \
-R DEGs.bed -b 5000 -a 5000 \
--outFileName matrix.1to4.gz \
--regionBodyLength 5000 --numberOfProcessors 5 --samplesLabel WT_1 WT_2 3a1KO_1 3a1KO_2

plotProfile -m ./deeptool_result/matrix.1to4.gz -out ./deeptool_result/density_DEGs.pdf --plotTitle "" \
--numPlotsPerRow 3 --perGroup --outFileNameData ./deeptool_result/density_DEGs_perGroup_5kb.tab
```

```
library(data.table)
avg.WGBS <- fread("density_DEGs_perGroup_5kb.tab",skip=1,fill=TRUE)
avg.WGBS <- avg.WGBS[-1,1:1502]
iOrd <- paste0(avg.WGBS$V1,"_",avg.WGBS$V2)
avg.WGBS <- as.data.frame(t(avg.WGBS[,3:1502]))
colnames(avg.WGBS) <- iOrd

dif.WGBS <- data.frame(UpDEGs=rowMeans(avg.WGBS[,c("3a1KO_1_UpDEGs","3a1KO_2_UpDEGs")]) - rowMeans(avg.WGBS[,c("WT_1_UpDEGs","WT_2_UpDEGs")]),
                       Other=rowMeans(avg.WGBS[,c("3a1KO_1_Other","3a1KO_2_Other")]) - rowMeans(avg.WGBS[,c("WT_1_Other","WT_2_Other")]),
                       DnDEGs=rowMeans(avg.WGBS[,c("3a1KO_1_DnDEGs","3a1KO_2_DnDEGs")]) - rowMeans(avg.WGBS[,c("WT_1_DnDEGs","WT_2_DnDEGs")])
)
rownames(dif.WGBS) <- 1:1500
df <- NULL
x <- data.frame(group="UpDEGs",value=dif.WGBS$UpDEGs,position=1:1500)
df <- rbind(df,x)
x <- data.frame(group="Other",value=dif.WGBS$Other,position=1:1500)
df <- rbind(df,x)
x <- data.frame(group="DnDEGs",value=dif.WGBS$DnDEGs,position=1:1500)
df <- rbind(df,x)

ggplot(df, aes(position, value, colour = group)) + 
  geom_point(size=0.1)+
  geom_smooth(method='loess', fullrange = TRUE,se=F,span=0.1) + 
  labs(y = "3a1KO-WT", x = "position")+
  theme_bw()+
  theme(legend.position = "top", legend.text = element_text(size = 12), 
        legend.title = element_text(size = 14)) + ylim(c(-0.125,0))
```

### Valcano plot of methylation changes for DMRs

```R
library(ggplot2);library(Seurat)
# load 
DMR_3a1KO <- fread("./MethyRatio/DMR/3a1KO_vs_WT.txt")
DMR_3a1KO <- data.frame(DMR_3a1KO,stringsAsFactors = F)
  df <- data.frame(WT=DMR_3a1KO$mean.g1,
                   V3a1KO=DMR_3a1KO$mean.g2,
                   pvalue=-1*log10(DMR_3a1KO$p..MWU.),
                   FDR = ifelse(DMR_3a1KO$q.value>=0.05,"","FDR<0.05"))
  p1 <- ggplot(df, aes(x=WT, y=V3a1KO, color=pvalue)) + geom_point(size=0.1) +
    scale_color_gradient2(midpoint=median(df$pvalue), low="#525252", mid="white",
                          high="red", space ="Lab" ) +
    theme_classic() 
  p2 <- AugmentPlot(p1, width = 10, height = 10, dpi = 300)
```



## Softwares

BSMAP (v2.90) 

metilene (v0.2-8) 

iBSTools (https://github.com/methylation/ iBSTools)

DeepTools (v3.1.3) 

R library “EnrichedHeatmap” (v1.14.0) 

Seurat (v4.0.3)

ggplot2 (v3.3.5)