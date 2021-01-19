This is the repository that contains the analysis of the lung adenocarcinoma single cell dataset

## Getting started

Clone the repo
Download the Data_input folder from the link below into the repo
https://drive.google.com/drive/folders/1nONsp9VuhmPzuDvMet0i8x26eV9r5lkT?usp=sharing 

## Representative ChIP-seq Scripts

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

### Using deeptools to normalize and visualize the data

```
module load deeptools
bamCoverage --bam $samplename --outFileName ./BWs/${samplename%.sort.bam}.bw -p 8 --effectiveGenomeSize 2150570000 --normalizeUsing RPGC --extendReads 200 --binSize 30 --outFileFormat bigwig

computeMatrix reference-point -S $sample.bw \
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


## Softwares

BWA-MEM algorithm (v0.7.17) 

SAMtools (v1.10) (http://www.htslib.org/). 

Phantompeakqualtools (https://github.com/kundajelab/phantompeakqualtools). 

DeepTools (v3.1.3) 