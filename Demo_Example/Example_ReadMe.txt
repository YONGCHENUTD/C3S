Here is an demo example.
The CAPTURE-3C-seq data is HS3 enhancer of K562 cell line that downloaded from NCBI GEO with access number of GSM2635075. The SRA number is SRR5583324. The central position of sgRNA target of this HS3 region is chr11:5305934.

### You can download it with following command under current folder.
> wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-7/SRR5583324/SRR5583324.1
### This needs to dump reads into two separated files for reads R1 and R2 (with SRA Toolkit installed https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc)
> fastq-dump --split-files SRR5583324.1 
### Optionally, a preprocessed data can be downloaded here 
Reads file1: SRR5583324_1.fastq.gz https://drive.google.com/file/d/1tHhPFmEZ0SUKHh4MoTK5vP3Xl5Mjp8G9/view?usp=sharing
Reads file2: SRR5583324_2.fastq.gz https://drive.google.com/file/d/1lN8ahayDZUzuO8auyzonIzoyacx7RdMF/view?usp=sharing

### You need download the genome version you need for bowtie2 mapping from bowtie2 website. For example, 
> wget ftp://ftp.ncbi.nlm.nih.gov/genomes/archive/old_genbank/Eukaryotes/vertebrates_mammals/Homo_sapiens/GRCh38/seqs_for_alignment_pipelines/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.bowtie_index.tar.gz 

### Now you have your data and genome available. You can run C3S command as following.
> runC3S.py -x hg38 -1 RR5583324_1.fastq.gz -2 RR5583324_1.fastq.gz sample_R2.fastq.gz --prefix HS3 --bait chr11:5305934

### After running, you will have results in current folder 
