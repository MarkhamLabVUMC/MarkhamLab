# MarkhamLab
A repo that contains exploratory data analysis based on a high-dimensional transcriptomic dataset.

## Contents
*Note it is recommended to go through these files in order.

### main Branch

**Weekly_Meeting_Notes.md:** Key points discussed during weekly meetings and the work to be done for the following week.

### kallisto_guide Branch

**Kallisto_Notes.md:** A high level overview of what the program is and the input it takes.

**Kallisto_Installation.md:** A step-by-step tutorial on how to install Kallisto on your device.

**Coding_Kallisto.ipynb:** A step-by-step tutorial on how to run Kallisto in Python and R, as well as descriptions of key commands.

### data_transformation Branch

**Data_Setup.md:** Guidelines on how to properly mount and eject the external hard drive in the WSL environment and overall system.

**Pseudoalign_Samples.ipynb:** A Jupyter Notebook that pseudoaligns the samples and maps them to a list of genes using Kallisto.

**Counts_Matrix_tximport.R:** A R file that produces the counts matrix for the samples.

**DESeq2_Analysis.md:** A step-by-step tutorial on how to do the DESeq2 tutorial.

**Volcano_Plots.R:** A R script that produces volcano plots with the DESeq pipeline.