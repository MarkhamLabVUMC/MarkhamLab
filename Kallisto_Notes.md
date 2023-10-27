# Kallisto Notes

This document explains the Kallisto library using multiple references.

**GitHub Page: (https://github.com/pachterlab/kallisto)**

- Program for quantifying abundances of transcripts from RNA-Seq data
    - More generally of target sequences using high-throughout sequencing reads
    - Based on idea of 'pseudoalignment'
        - Rapidly determining compatibility of reads with targets, without need for alignment
    - Quantify 30 million human bulk RNA-Seq reads in less than 3 minutes on Mac desktop computer
        - Using only read sequences and transcriptome index that takes 10 minutes to build
    - Pseudoalignment of reads preserves key information needed for quantification

**ChatGPT**

Prompt: I am a data science grad student. I want to learn about the kallisto library and what it does. Please explain it to me nicely and simply as I do not know anything about it.

- Program for quantifying abundances of transcripts from RNA-Seq data
- RNA-Seq is technique used to profile expression of genes in biological sample
    - Lets you see which genes are "turned on" or "turned off" in a sample, and to what degree
- Quantifying abundances of transcripts means counting number of transcripts for each gene
    - Each gene can produce a molecule called transcript (mRNA)
    - Number of transcripts can tell us about activity of gene
        - More transcripts indicates gene is more active
- Kallisto uses 'pseudoalignment'
    - Traditional way is to align reads directly to reference genome
        - Too time-consuming and intense to compute
    - Kallisto maps them to set of reference transcripts
        - Method is faster and less intense to compute
        - Idea is to determine compatibility of reads with transcripts, without needed full alignment
- Important because can allow researchers to quickly and accurately estimate transcript abundances from RNA-Seq data
    - Speed up research projects, gain insight faster
- Key features:
    - Speed
        - Process millions of reads in minutes
    - Accuracy
        - Produced results are comparable in accuracy to other methods that take quite longer
    - Bootstrap Sampling
        - Has built-in feature for bootstrap resampling
        - Provides measures of uncertainty is estimated abundances

Prompt: Can you please elaborate on the functionality of kallisto? As in, certain functions and a beginner level tutorial on how to use it?

- Key Functionalities:
    - Indexing
        - Before quantifying transcripts, need to build index from set of reference transcripts
            - 1 time process for each set
    - Quantification
        - Once index is built, can quantify transcript abundances for any RNA-Seq dataset using that index
    - Pseudoalignment
        - Underlying method that allows Kallisto to be quick
        - Instead of aligning each read to reference, determines compatibility of read with each transcript
    - Bootstrap Sampling
        - Can perform bootstrap sampling to provide estimates of uncertainty associated with each transcript's abundance
- Tutorial
    1. Installation
        - Need to download from Kallisto GitHub repository
        - Follow provided installation instructions
    2. Indexing
        - Need to have a set of reference transcripts in FASTA format to build index
        - Ex: kallisto index -i transcripts.idx transcripts.fasta
            - '-i' specifies name of output index file
            - 'transcripts.idx' is name of output index file
            - 'transcripts.fasta' is input FASTA file containing reference transcripts
    3. Quantification
        - Once index is built, can quantify transcript abundances from RNA-Seq data
            - Ex: Assume you have paired-end RNA-Seq data in 2 files 'reads_1.fastq' and 'reads_2.fastq'
                - kallisto quant -i transcripts.idx -o output_folder reads_1.fastq reads_2.fastq
                    - '-i' specifies name of index file
                    - 'o' specifies name of output folder where results will be saved
                    - 'reads_1.fastq' and 'reads_2.fastq' is input RNA-Seq data files
    4. Results
        - After 'quant' command is run. output is several files in 'output_folder'
        - Primary file of interest us 'abundance.tsv'
            - Contains estimated counts and TPM (Transcripts Per Million) values for each transcript
- Helpful Tips
    - Always ensure reference transcripts match organism and version used in your experiment
    - If unsure about any command, use 'kallisto help' to see list of available commands and options
    - Read Kallisto documentation to gain deeper understanding


