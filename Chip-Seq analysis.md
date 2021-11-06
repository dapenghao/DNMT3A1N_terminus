Representative ChIP-seq Scripts

### Mapping of ChIPSeq data using BWA

```bash
module load bwa
gunzip -cd $samplename > ${samplename%.gz}
samplename=${samplename%.gz}
bwa mem -t 8 -c 1 ./mm9/mm9.fa $samplename > ${samplename%.fastq.gz}.sam
```




### Generating bam files using SAMtools

```bash
module load samtools

# convert Sam files to sorted bam files:
samtools view --threads 8 -bS $samplename > ../bams/${samplename%.sam}_1.bam
# remove duplicates
cd ../bams
samtools rmdup -s ${samplename%.sam}_1.bam ${samplename%.sam}.rmdup.bam
# sort bams
samtools sort --threads 8 ${samplename%.sam}.rmdup.bam -o ${samplename%.sam}.sort.bam
# index bams
samtools index ${samplename%.sam}.sort.bam
```

### Using phantompeakqualtools to assess enrichment and fragment length

```bash
Rscript ./Tools/phantompeakqualtools/run_spp.R -c=sample.sort.bam -savp=./${i%.sort.bam}.pdf -out=./xcor/xcor_metrics_${i%.sort.bam}.txt
```

### Using deeptools to normalize data

```bash
module load deeptools
bamCoverage --bam $samplename --outFileName ./BWs/${samplename%.sort.bam}.bw -p 8 --effectiveGenomeSize 2150570000 --normalizeUsing RPGC --extendReads 200 --binSize 30 --outFileFormat bigwig
```

### Input subtraction

```bash
module load deeptools
bigwigCompare -b1 3a1KO.bw -b2 input.bw \
--binSize 30 \
--outFileFormat bigwig --outFileName 3a1KO_inputSubtract.bw \
--operation subtract -p 8
```

### Data visualization

```
module load deeptools
computeMatrix reference-point -S 3a1KO_inputSubtract.bw \
-R Protein-coding-genes.bed \
-b 5000 -a 5000 --outFileName ./matrix.tss_sample.gz \
--referencePoint=TSS --numberOfProcessors 5 --binSize 50

plotHeatmap --matrixFile ./matrix.tss_sample.gz \
--outFileName ./TSS_sample.pdf --colorMap Reds \
--zMin 0 --zMax 5 --whatToShow 'plot, heatmap and colorbar' --sortRegions keep

plotProfile --matrixFile matrix.tss_sample.gz \
--outFileName ./signal_ES_sample.pdf \
--yMin 0.3 --yMax 2.2 --numPlotsPerRow 1 --perGroup
```

### Peak calling

```bash
module load macs2/2.1.2
macs2 callpeak -f BAM -t H3K4me3.sort.bam -c input.sort.bam -n H3K4me3 --keep-dup all -g mm --broad --broad-cutoff 0.05 --nomodel --extsize 200 --max-gap 500 --min-length 500
```

### Genome-wide correlations

```bash
module load deeptools
multiBigwigSummary bins \
-b sample_1.bw \
sample_2.bw \
sample_3.bw \
-out scores_per_bin_ES.npz

plotCorrelation -in scores_per_bin_ES.npz \
--corMethod spearman --skipZeros \
--plotTitle "spearman Correlation of Genome-wide 10kb bins" \
--whatToPlot heatmap --colorMap RdBu_r --plotNumbers \
-o heatmap_SpearmanCorr_ES.pdf \
--outFileCorMatrix SpearmanCorr_ES.tab --removeOutliers \
--blackListFileName mm9-blacklist.bed
```

### Identification of enhancers

```bash &amp; R
module load bedtools
#1.report those entries in H3K4me1 that have _no overlaps_ with H3K4me3:
bedtools intersect -a H3K4me1_peaks.bed -b H3K4me3_peaks.bed -v > Enhancers.bed

#2.download TSS of mm9 known genes from https://ccg.epfl.ch/mga/mm9/ucsc/ucsc.html and convert to bed file:
#Table was downloaded with the following attributes: mm9.knownGene, mm9.ensemblSource, mm9.kgXref, mm9.knownToEnsembl, mm9.refSeqStatus, mm9.spMrna.
#Genes were retained if annotated as "protein coding" or if they had a protein ID or a refSeq status.
# TSS_1kb_mm9.bed was generated for 1kb flanking region of genes' TSS.

#3.report those entries in Enhancers that have _no overlaps_ with TSS ± 1 kb:
bedtools intersect -a Enhancers.bed -b TSS_mm9.bed -v > Enhancers2.bed

#4. Neighboring enhancers located within 500 bp were merged.
bedtools merge -i Enhancers2.bed -d 500 > Enhancers.bed

#5. Enhancers were categorized into active and poised  according to overlapping with H3K27ac peaks or not, respectively:
#Perform a "left outer join"
bedtools intersect -a Enhancers.bed -b H3K27ac_peaks.bed -loj > Enhancers2.bed
# according the output, define active and poised enhancers
```



## Softwares

BWA-MEM algorithm (v0.7.17) 

SAMtools (v1.10) (http://www.htslib.org/). 

Phantompeakqualtools (https://github.com/kundajelab/phantompeakqualtools). 

DeepTools (v3.1.3) 

MACS2 (v2.1.2)

Bedtools (v2.28.0)