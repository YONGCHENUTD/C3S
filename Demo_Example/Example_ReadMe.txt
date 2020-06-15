Here is a demo example. The CAPTURE-3C-seq data of HS3 enhancer of K562 cell line was downloaded from NCBI GEO with access number of GSM2635075. The SRA number is SRR5583324. The central position of sgRNA target of this HS3 region is chr11:5305934.

You can download the SRA file using prefetch. Then dump the fastq reads from the SRA file. (SRA Toolkit: https://github.com/ncbi/sra-tools).
> prefetch SRR5583324
> fastq-dump --split-3 --helicos --gzip  ~/ncbi/public/sra/SRR5583324.sra

In case of difficulty in installation of the SRA Toolkit, the fastq reads can be downloaded from Googl Drive:
> wget "https://drive.google.com/file/d/1tHhPFmEZ0SUKHh4MoTK5vP3Xl5Mjp8G9/view?usp=sharing" -O SRR5583324_1.fastq.gz
> wget "https://drive.google.com/file/d/1lN8ahayDZUzuO8auyzonIzoyacx7RdMF/view?usp=sharing" -O SRR5583324_2.fastq.gz

Download the human genome sequence file, and get rid of the non-refernece chromosomes. Then build the bowtie2 index.
> wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
> gunzip hg19.fa.gz
> samtools faidx hg19.fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY >hg19_basic.fa
> bowtie2-build hg19_basic.fa hg19_basic

Now you can run C3S using the folloing command:
> runC3S.py -x hg19_basic -1 SRR5583324_1.fastq.gz -2 SRR5583324_2.fastq.gz --prefix HS3 --bait chr11:5305934

Result files:
/010ReadMapping:
  ├── HS3_R1_bowtie2.log
  ├── HS3_R1_un.fastq.gz
  ├── HS3_R1_samtools.log
  ├── HS3_R1.bam
  ├── HS3_R1_flagstat.log
  ├── HS3_R2_bowtie2.log
  ├── HS3_R2_un.fastq.gz
  ├── HS3_R2_samtools.log
  ├── HS3_R2.bam
  ├── HS3_R2_flagstat.log
  ├── HS3_R1_split.fastq.gz
  ├── HS3_R2_split.fastq.gz
  ├── HS3_R1_remap_bowtie2.log
  ├── HS3_R1_remap_un.fastq.gz
  ├── HS3_R1_remap_samtools.log
  ├── HS3_R1_remap.bam
  ├── HS3_R1_remap_flagstat.log
  ├── HS3_R2_remap_bowtie2.log
  ├── HS3_R2_remap_un.fastq.gz
  ├── HS3_R2_remap_samtools.log
  ├── HS3_R2_remap.bam
  ├── HS3_R2_remap_flagstat.log
  ├── HS3.pairs.gz
  ├── HS3.pairs.gz.tbi
/020Plotting:
  ├── HS3_stats.pdf
/030Model:
  ├── HS3_wu.bedpairs
  
