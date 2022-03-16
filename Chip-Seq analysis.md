## 1. Processing of ChIP-Seq data

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



## 2. Data visualization

### Heatmap of Signals around TSS

Related to "Figure.2a, 6a, 7a, 7b, 7h; Extended Data Figure.9d"

```bash
module load deeptools
computeMatrix reference-point -S S1_inputSubtract.bw S2_inputSubtract.bw\
-R Protein-coding-genes.bed \
-b 5000 -a 5000 --outFileName ./matrix.tss_sample.gz \
--referencePoint=TSS --numberOfProcessors 5 --binSize 50

plotHeatmap --matrixFile ./matrix.tss_sample.gz \
--outFileName ./TSS_sample.pdf --colorMap Reds \
--zMin 0 --zMax 5 --whatToShow 'heatmap and colorbar' --sortRegions keep
```

### Heatmap of Signals around Gene Body

Related to "Figure.5f, 5h; Extended Data Figure.3a, 3i, 8h"

```bash
module load deeptools
computeMatrix scale-regions -S S1_inputSubtract.bw S2_inputSubtract.bw \
-R Protein-coding-genes.bed \
-b 2000 -a 2000 --outFileName ./matrix.Genes_sample.gz \
--regionBodyLength 5000 --numberOfProcessors 7 --binSize 50 --missingDataAsZero

plotHeatmap --matrixFile ./matrix.Genes_sample.gz \
--outFileName ./Genes_sample.pdf --colorMap Reds \
--zMin 0 --zMax 5 --whatToShow 'heatmap and colorbar' --sortRegions keep
```

### Heatmap of Signals around defined locations

#### Heatmap of Signals around H2AK119ub Peaks

Related to "Figure.7i"

###### Peak calling

```bash
module load macs2/2.1.2
macs2 callpeak -f BAM -t ./H2AK119ub.sort.bam -c ./input.sort.bam -n H2AK119ub_ES --keep-dup all -g mm --broad --broad-cutoff 0.05 --nomodel --extsize 135 --max-gap 500 --min-length 500
```

###### Heatmap of signals

```bash
computeMatrix reference-point -S \
WTJ1_H2AK119ub_inputSubtract.bw \
S1_inputSubtract.bw S2_inputSubtract.bw ... Sn_inputSubtract.bw\
-R ./MACS/H2AK119ub_ES_peaks.bed \
-b 10000 -a 10000 --outFileName ./matrix.peaks_ES.gz \
--referencePoint=center --numberOfProcessors 5 --binSize 50 \
--missingDataColor 0 --sortUsingSamples 7 --sortRegions descend

plotHeatmap --matrixFile ./matrix.peaks_ES.gz \
--outFileName ./H2AK119ub_Peaks_ES_inputSubtracted_10kbFlanking.pdf --colorMap Reds \
--zMin 0 --zMax 3 --whatToShow 'heatmap and colorbar' --sortRegions keep
```

#### Signals around Enhancers

Related to "Figure.; Extended Data Figure.8e"

###### Identification of enhancers

```bash
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

#5. Enhancers were categorized into active and poised according to overlapping with H3K27ac peaks or not, respectively:
#Perform a "left outer join"
bedtools intersect -a Enhancers.bed -b H3K27ac_peaks.bed -loj > Enhancers2.bed
# according the output, define active and poised enhancers
```

###### Heatmap around enhancers

```bash
module load deeptools
computeMatrix reference-point -S \
3aFLAG_H3K4me1_inputSubtract.bw \
3aFLAG_H3K27ac_inputSubtract.bw \
3aFLAG_H3K4me3_inputSubtract.bw \
3aFLAG_H3K27me3_inputSubtract.bw \
3A1.bw \
-R Enhancers2.bed \
-b 5000 -a 5000 --outFileName ./matrix.Enhancers_Cortex.gz \
--referencePoint=center --numberOfProcessors 5 --binSize 50 \
--missingDataAsZero --sortUsingSamples 3 --sortRegions descend

plotHeatmap --matrixFile ./matrix.Enhancers_Cortex.gz \
--outFileName ./Enhancer_inputSubtracted_byK4me1.pdf --colorMap Reds \
--zMin 0 --zMax 5 --whatToShow 'heatmap and colorbar' --sortRegions keep \
--heatmapWidth 5
```

### Density Plot of Signals around TSS or Gene Body

Related to "Figure.2b, 6c, 7c, ; Extended Data Figure.3i, 9c"

```bash
module load deeptools
#1. Around TSS:
computeMatrix scale-regions -S S1_inputSubtract.bw S2_inputSubtract.bw \
-R Protein-coding-genes.bed \
-b 2000 -a 2000 --outFileName ./matrix.tss_sample.gz \
--regionBodyLength 5000 --numberOfProcessors 7 --binSize 50 --missingDataAsZero

plotProfile --matrixFile matrix.tss_sample.gz \
--outFileName ./signal_tss_sample.pdf \
--yMin -1 --yMax 1 --numPlotsPerRow 1 --perGroup
#2. Around Gene Body:
computeMatrix scale-regions -S S1_inputSubtract.bw S2_inputSubtract.bw \
-R Protein-coding-genes.bed \
-b 2000 -a 2000 --outFileName ./matrix.Genes_sample.gz \
--regionBodyLength 5000 --numberOfProcessors 7 --binSize 50 --missingDataAsZero

plotProfile --matrixFile matrix.Genes_sample.gz \
--outFileName ./signal_Genes_sample.pdf \
--yMin -1 --yMax 1 --numPlotsPerRow 1 --perGroup
```



## 3. Genome-wide correlations

Related to "Extended Data Figure.9i"

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



## 4. Softwares

BWA-MEM algorithm (v0.7.17) 

SAMtools (v1.10) (http://www.htslib.org/). 

Phantompeakqualtools (https://github.com/kundajelab/phantompeakqualtools). 

DeepTools (v3.1.3) 

MACS2 (v2.1.2)

Bedtools (v2.28.0)

