Here is an demo example.
The CAPTURE-3C-seq data is HS3 enhancer of K562 cell line that downloaded from NCBI GEO with access number of GSM2635075. The SRA number is SRR5583324. The central position of sgRNA target of this HS3 region is chr11:5305934.

### You can download it with following command under current folder.
> wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-7/SRR5583324/SRR5583324.1
### This needs to dump reads into two separated files for reads R1 and R2 (with SRA Toolkit installed https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc)
> fastq-dump --split-files SRR5583324.1 
### Optionally, a preprocessed data can be downloaded here 
Reads file1: SRR5583324_1.fastq.gz https://drive.google.com/file/d/1tHhPFmEZ0SUKHh4MoTK5vP3Xl5Mjp8G9/view?usp=sharing
Reads file2: SRR5583324_2.fastq.gz https://drive.google.com/file/d/1lN8ahayDZUzuO8auyzonIzoyacx7RdMF/view?usp=sharing

### You need download the genome version you need for bowtie2 mapping from bowtie2 website. For example, to build GRCh38.
### download the genome sequence file:
> wget ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_34/GRCh38.primary_assembly.genome.fa.gz
> gunzip GRCh38.primary_assembly.genome.fa.gz

### Extract the reference genomes
> samtools faidx GRCh38.primary_assembly.genome.fa chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 chr20 chr21 chr22 chrX chrY >hg38_ref.fa

### Build bowtie2 index
> bowtie2-build hg38_ref.fa hg38_ref

### Now you have your data and genome available. You can run C3S command as following.
> runC3S.py -x hg38_ref -1 SRR5583324_1.fastq.gz -2 SRR5583324_1.fastq.gz --prefix HS3 --bait chr11:5305934

### After running, you will have results in current folder 
/010ReadMapping:
HS3_R1_bowtie2.log
HS3_R1_un.fastq.gz
HS3_R1_samtools.log
HS3_R1.bam
HS3_R1_flagstat.log
HS3_R2_bowtie2.log
HS3_R2_un.fastq.gz
HS3_R2_samtools.log
HS3_R2.bam
HS3_R2_flagstat.log
HS3_R1_split.fastq.gz
HS3_R2_split.fastq.gz
HS3_R1_remap_bowtie2.log
HS3_R1_remap_un.fastq.gz
HS3_R1_remap_samtools.log
HS3_R1_remap.bam
HS3_R1_remap_flagstat.log
HS3_R2_remap_bowtie2.log
HS3_R2_remap_un.fastq.gz
HS3_R2_remap_samtools.log
HS3_R2_remap.bam
HS3_R2_remap_flagstat.log
HS3.pairs.gz
HS3.pairs.gz.tbi

/020Plotting:
HS3_stats.pdf

/030Model:
HS3_wu.bedpairs


