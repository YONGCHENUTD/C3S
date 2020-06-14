**C3S: model-based analysis and pipeline of dCas9 Capture-3C-Seq data**

C3S is a model-based analysis and pipeline of dCas9 Capture-3C-Seq data (Xin Liu et.al Cell 2017). The CAPTURE-3C-seq method now is widely used and highly contributed to studies of the cancer epigenetics and other fundamental biological questions. Here C3S is introduced to analyse CAPTURE-3C-seq data from raw sequencing reads to significantly interacted chromatin loci. It uses multiplescale Bayesian models for the significance calling of chromatin interactions, especialluy providing different models to analyse intra- and inter-chromosomal chromatin interactions for different mammalian species. It oprovides a versatile and flexible pipeline to analyze the dCas9 Capture-3C-Seq data.

=============================

1. Prerequisition
-------------------
- Python 2.7 packages (automatically installed)

  - numpy >= 1.13.1
  - scipy >= 0.19.1
  - statsmodels >=0.8.0
  - pandas >= 0.15.2
  - matplotlib >= 2.0.2
  - seaborn >= 0.7.1
  - pysam >= 0.11.2.2

- Other tools

  - bowtie2 >= 2.2.2
  - samtools >= 1.5
  
2. Installation
----------------

::

  > git clone https://github.com/YONGCHENUTD/C3S.git
  > cd C3S
  > python setup.py install --user

3. Usage of C3S
----------------

- Usage of runC3S.py

::

  usage: runC3S.py -x hg38 -1 sample_R1.fastq.gz [sample_R1.fastq.gz ...] -2
                 sample_R2.fastq.gz [sample_R2.fastq.gz ...] --prefix prefix
                 [--bait chr11:5305934] [--extendsize 100000] [--readlen 36]
                 [--seed 1024] [--smooth-window 100] [--nperm 10000] [-w "."]
                 [-p 10]

- Required parameters:

+--------------------------------------+--------------------------------------------------------------+
|-x hg38, --genome hg38                |Bowtie2 built genome.                                         |
+--------------------------------------+--------------------------------------------------------------+
|-1 S1_R1.fastq.gz [S2_R1.fastq.gz ...]|Read 1 fastq file. Can be gzip(.gz) or bzip2(.bz2) compressed.|
+--------------------------------------+--------------------------------------------------------------+
|-1 S1_R2.fastq.gz [S2_R2.fastq.gz ...]|Read 2 fastq file. Can be gzip(.gz) or bzip2(.bz2) compressed.|
+--------------------------------------+--------------------------------------------------------------+
|--prefix                              |Prefix of result files.                                       |
+--------------------------------------+--------------------------------------------------------------+
|--bait chr position                   |Bait genomic locus. [For example "chr11:5305934"]             |
+--------------------------------------+--------------------------------------------------------------+

- Optional parameters:

+--------------------------------------+--------------------------------------------------------------+
|--extendsize 100000                   |Length to be extended from bait regions. [Defaut=100000].     |
+--------------------------------------+--------------------------------------------------------------+
|--readlen 36                          |Read length. [Default=36]                                     |
+--------------------------------------+--------------------------------------------------------------+
|--seed 1024                           |Seed to generate random values. [Default=1024].               |
+--------------------------------------+--------------------------------------------------------------+
|-smooth-window 101                    |Smooth window for peak size inference. [Default=101].         |
+--------------------------------------+--------------------------------------------------------------+
|--nperm 10000                         |Number of permutatons. [Default=10000].                       |
+--------------------------------------+--------------------------------------------------------------+
|-w "."                                |Working directory. [Default="."].                             |
+--------------------------------------+--------------------------------------------------------------+
|-p 10                                 |Number of processes. [Default=10].                            |
+--------------------------------------+--------------------------------------------------------------+


4. Examples:
-----------------

::

Here is an demo example. The CAPTURE-3C-seq data is HS3 enhancer of K562 cell line that downloaded from NCBI GEO with access number of GSM2635075. The SRA number is SRR5583324. The central position of sgRNA target of this HS3 region is chr11:5305934.

- Download the SRA file using prefetch. The downloaded file is in ~/ncbi/public/sra by default. (SRA Toolkit: https://github.com/ncbi/sra-tools). 

::

> prefetch SRR5583324

- Dump fastq reads from SRR5583324.sra

::

> fastq-dump --split-3 --helicos --gzip  ~/ncbi/public/sra/SRR5583324.sra


- In case of difficulty in installation of the SRA Toolkit, the fastq reads can be downloaded from Googl Drive: 

::

> wget "https://drive.google.com/file/d/1tHhPFmEZ0SUKHh4MoTK5vP3Xl5Mjp8G9/view?usp=sharing" -O SRR5583324_1.fastq.gz
> wget "https://drive.google.com/file/d/1lN8ahayDZUzuO8auyzonIzoyacx7RdMF/view?usp=sharing" -O SRR5583324_2.fastq.gz


- Download the human genome sequence file, and get rid of the non-refernece chromosomes.

::

> wget http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
> gunzip hg19.fa.gz
> samtools faidx hg19.fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY >hg19_basic.fa

- Build bowtie2 index

::

> bowtie2-build hg19_basic.fa hg19_basic


- Run C3S using the folloing command:

::

> runC3S.py -x hg19_basic -1 SRR5583324_1.fastq.gz -2 SRR5583324_2.fastq.gz --prefix HS3 --bait chr11:5305934


- Result files: 

::

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
    |-HS3_stats.pdf
  /030Model:
    |-HS3_wu.bedpairs

5. Citations of C3S
----------------------------------

::

Yong Chen, Yunfei Wang, Xin Liu, Jian Xu, Michael Q. Zhang. Model-based Analysis of Chromatin Interactions from dCas9-Based CAPTURE-3C-seq. PLOS ONE

Liu X, Zhang Y, Chen Y, et al. In Situ Capture of Chromatin Interactions by Biotinylated dCas9. Cell. 2017;170(5):1028‐1043.e19. doi:10.1016/j.cell.2017.08.003


  
