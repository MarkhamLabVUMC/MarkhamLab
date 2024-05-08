# Weekly Meeting Notes

**October 20, 2023**

- Colon cancer exists in middle of pipe
    - Loaded with bacteria
- Some colon tissue & colon cancer tissue covered in bacteria
- Big question: Relationship between bacteria touching colon cancer & colon cancer itself?
- Looking carefully at cancerous tissue
    - To see if there is any bacteria
- Epithelial cells: Cells those lines tissue
    - For any organ that secrets things
        - Breast, prostate
    - Most cancers come from epithelial cells
    - Separates us from microbes inside of us
        - Special, thin barrier
            - Ability to absorb nutrients, excrete waste, keep out pathogens
- In colon cancer & in some normal tissue, some bacteria invade or get very close to epithelial cells
    - Biofilm
        - Sugary protein that is goo that helps bacteria adhere to surface
        - Help protect bacteria from antibiotics
        - Resistant to own immune system
        - More common in right side of colon
        - Ask question: Which came first? Bacteria nestle up next to epithelial cells to form cancer or cancer formed first and bacteria got attracted
- C. difficile was a part of these biofilms
    - Bacteria that has some properties
        - Forms spores released in stool, carried onto next host (ex: via surfaces)
        - Spores in healthy GI tract do not do anything
            - In sick patients, will vegetate and produce toxins
                - TcdA, TcdB
                    - Mediate disease, kill epithelial cells
    - Patients who have C. difficile infection
        - Cause is these toxins
        - Become susceptible while taking antibiotics for something else
            - Can go into gut and kill healthy bacteria
    - Was not previously identified was related to causing colon cancer
        - Can accelerate colon cancer in genetically susceptible mice
- Question: How does C. difficile do that?
    - If toxin is what mediate disease in infectious diarrhea, is toxin what in a different patient lead to cancer?
        - Wanted to test that where can remove several variables
            -How big is person, diet, bathroom frequency
- Test it out in dishes
    - Matrigel
    - Expose them just to toxin, taking out many variables from C. difficile
        - Is it possible for C. difficile toxin itself to effect epithelial cells to induce cancer?
- 3728T
    - Identification of sample patient, T for tumor
- Consortium from 30 bacteria derived from original slurry from patient 3728 or mice who got slurry 3728
- At step where have raw sequencing data
    - Each of cDNA got sequenced in both directions
    - During reverse transcription single strand mRNA becomes into double strand DNA
    - Parameters built into it
- 2 separate files for samples
    - Samples numbered 1-30
    - Each have R1 and R2 file
        - Annotated by project ID, initials of person, sample number, R1 & R2 for separate files (fastQ files, compressed with gzip)
- Storage & Memory
    - Memory
    - Align these 1 by 1, computationally intense, require way more memory
    - Use package: kallisto 
- Next Steps: 
    - Look through documentation of kallisto
    - Need to set up private GitHub
        - Find out how much data can be uploaded

**October 27, 2023**
- Showed Dr. Markham GitHub Repo and tutorials created for kallisto
- Next Steps: 
    - Look at 'sleuth' package in R
    - Move raw data onto external hard drive for access
    - Experiment importing & reading in data to Kallisto
        - Follow manual to do alignments & mapping

**November 3, 2023**
- No meeting as data is currently being transferred to external hard drive

**November 10, 2023**
- Data
    - fastq file format
    - R1 and R2, 2 directions for each sample
    - Compressed, hence fastq.gz
    - Used illumina sequencing technology
        - Paired-end run
        - Novaseq 6000 machine was used to do sequencing
    - Another 15 files
- Summary excel sheet gives more information
    - Total Yield
    - Q30: Measure of quality
- Next Steps:
    - Look into what FASTQ files are
    - Try to pseudoalign 1 sample & map them to list of genes and list of counts for each gene
        - R1 & R2
    - Having genes as rows, columns representing samples 1-30

**November 17, 2023**
- No meeting, continue working on pseudoalignment

**November 24, 2023**
- No meeting, Thanksgiving break

**December 1, 2023**
- Use est_counts, not tpm
- Normalize, transform data
- Get rest of samples on hard drive
- Next Steps:
    - After getting dataframe, take them & make comparisons across different samples
    - Group some samples together based on how they were treated experimentally
    - Statistical testing to compare genes across these groups
    - 0's in data is significant, leave it as is (Sparse matrices)
    - Programs that we will work with can handle these 0's (scanpy, based off of numpy and pandas)
    - Dataframe needs to be filtered, as some genes are 0 for some samples but not 0 for others, and then some samples that have 0 for all of them
    - Python packages for filtering RNA-Seq data

**Winter Break**
- Use new indexing
- Continue project with new index 

**January 8 - February 5, 2024**
- Quantified samples 1 - 18, 20 - 30 with new index

**February 7, 2024**
- Troubleshot quantification of sample 19

**February 9, 2024**
- Ensured all samples produced correct output
- Next Steps:
    - Need to start analysis of differently expressed genes (go from kallisto to counts matrix)
    - Generate counts matrix (Trinity via R or find a way in Python)
    - Produce a scatterplot of log2 fold changes against the mean normalized counts for all genes
    - Use PyDESeq2 for statistical analysis
    - Take statistically analyzed data and generate plots (using packages like seaborn or matplotlib)

**February 15, 2024**
- Discussed Trinity was not working to generate counts matrix
- Showed alternate counts matrix, was not desired output
- Next Steps:
    - Delete every other target_id column from alternate counts matrix as they are identical (only need 1 in total)
    - Need to change target_id to names of genes using names from ensemble database
    - Try tximport package for counts matrix

**February 23, 2024**
- Showed new counts matrix that was a result of tximport (in R)
- Made minor tweaks to get final counts matrix in .csv
- Next Steps:
    - Use counts matrix as input to finish rpub tutorial
    - Switch to R for analysis and do DESeq2 tutorial

**March 1, 2024**
- Working on DESeq2 tutorial
- Next Steps:
    - Finish DESeq2 tutorial

**March 22, 2024**
- Showed progress on DESeq2 tutorial
    - Heatmaps
- Next Steps:
    - Plot a histogram, where the x-axis is the counts (0 - 1000), y-axis is # of genes
        - Plot all conditions in 1 plot
    - Figure out how to get row names visible on heatmap
    - Use 'Condition' and 'Cell' when comparing 2 variables in DESeq2 analysis
    - Create PCA plot

**March 29, 2024**
- Accomplished all previous tasks that were requested
- Made more progress on DESeq2 analysis
- Clarified what parts of the analysis are not necessary
- Next Steps:
    - Create a heatmap of the 20 most differentially (regulated) expressed genees
    - Modify the PCA plot so that the shape of the points vary for the Cell type
    - Produce presentation worthy (informative) plots for Dr. Markham's talk

**April 5, 2024**
- Discussed possibility of continuing project next semester as a capstone
    - Showed document that needs to be filled out
- Next Steps:
    - Fill out capstone document
    - Continue to work on DESeq2 analysis

**April 12, 2024**
- Finished coding portion of DESeq2 analysis
- Showed capstone document
- Next Steps:
    - Render DESeq2 analysis into .rmd so it is easier to follow along
    - Try creating volcano plots
    - Send capstone document for dean's approval

**April 19, 2024**
- Wrapped up .rmd of DESeq2 analysis
- Showed various volcano plots produced by several different packages
- Next Steps:
    - Push everything to GitHub
    - Think about how to start capstone project

**April 23, 2024**
- Handed over external hard drive

**Summer Break**
- Cleaned up GitHub so summer intern can easily follow along the created pipeline
- Brainstormed next steps for capstone
    - Rough implementation, including what tools to use