Pipeline for global var gene expression estimation

Main pipeline takes paired-end RNA-seq (fq or fa) 

Pipeline is run in three steps:
1. Identify reads containing the LARSFADIG motif
2. De novo assemble the LARSFADIG (performed using rnaSPAdes)
3. Determine coverage over the middle of the LARSFADIG motif



This is executed using several multilanguage scripts and third party tools. Before you start, please ensure you have Python, seqtk, bwa tools, blastn, rnaSPAdes and samtools installed
