Here is an demo example.
The CAPTURE-3C-seq data is HS3 enhancer of K562 cell line that downloaded from NCBI GEO with access number of GSM2635075. The SRA number is SRR5583324.

### You can download it with following command
> wget https://sra-downloadb.be-md.ncbi.nlm.nih.gov/sos2/sra-pub-run-7/SRR5583324/SRR5583324.1
### This is to dump reads into two separated files for reads R1 and R2 (with SRA Toolkit installed https://trace.ncbi.nlm.nih.gov/Traces/sra/sra.cgi?view=toolkit_doc)
> fastq-dump --split-files SRR5583324.1 
### A preprocessed data can be downloaded here 
Reads file1: SRR5583324_1.fastq.gz https://drive.google.com/file/d/1tHhPFmEZ0SUKHh4MoTK5vP3Xl5Mjp8G9/view?usp=sharing
Reads file2: SRR5583324_2.fastq.gz https://drive.google.com/file/d/1lN8ahayDZUzuO8auyzonIzoyacx7RdMF/view?usp=sharing


### Now you can run C3S as following.
> 


