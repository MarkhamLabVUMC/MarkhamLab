# 🌟 MarkhamLab Repository
Welcome to the **MarkhamLab** repository! This project involves the comprehensive analysis of a high-dimensional transcriptomic dataset. From pseudoaligning RNA-Seq reads to building an interactive dashboard with gene-level insights, this repository documents every step of the process.

---

## 📂 Contents
Here’s a breakdown of the key components in this repository. **It is recommended to go through these files in order.**

### 🗓 **Weekly_Meeting_Notes.md**
- Summarizes the key points discussed during weekly meetings.
- Outlines tasks and deliverables for the following week, keeping the project on track.

### 🧬 **Kallisto_Notes.md**
- Provides a high-level overview of Kallisto, including its purpose and input requirements.
- Useful for beginners to understand how pseudoalignment works.

### 🛠 **Kallisto_Installation.md**
- Step-by-step instructions for installing Kallisto on your local device.
- Ensures a smooth setup for first-time users.

### 📓 **Coding_Kallisto.ipynb**
- A Jupyter Notebook showcasing how to run Kallisto using Python and R.
- Includes detailed explanations of key commands and use cases.

### 💾 **Data_Setup.md**
- Instructions on properly mounting and ejecting the external hard drive within the WSL environment.
- Covers system-specific guidelines to handle dataset storage.

### 🎯 **Pseudoalign_Samples.ipynb**
- A Jupyter Notebook demonstrating how to pseudoalign RNA-Seq samples and map them to a list of genes using Kallisto.
- Produces the data necessary for downstream analyses.

### 📊 **Counts_Matrix_tximport.R**
- An R script that processes pseudoaligned data to create a counts matrix for the samples.
- Uses the `tximport` package for importing and summarizing transcript-level data.

### 🔬 **DESeq2_Analysis.md**
- A comprehensive tutorial on running the DESeq2 pipeline.
- Explains every step, from pre-filtering data to performing differential expression analysis.

### 🌋 **Volcano_Plots.R**
- An R script for generating volcano plots based on DESeq2 results.
- Highlights significant genes, helping identify biological insights.

### 🔍 **Gene_Scraping.ipynb**
- A Python notebook for web scraping summaries and descriptions of human genes from the NCBI database.
- Prepares gene metadata for integration into downstream visualizations.

### 📊 **app.R**
- An R script powering the **Markham Lab DESeq Dashboard**.
- Features:
  - An **interactive chatbot** that allows users to query gene summaries and descriptions.
  - **Dynamic visualizations**, including volcano plots, MA plots, PCA plots, and heatmaps, for exploring transcriptomic data.
  - User-friendly filters to explore genes and conditions.

---

## 🛠 Technologies Used
- **Kallisto**: For pseudoalignment and transcript quantification.
- **R and RStudio**: For statistical analysis (DESeq2) and data visualization.
- **Python and Jupyter Notebooks**: For data preprocessing, scraping, and exploratory analysis.
- **Shiny**: To create an interactive dashboard for gene-level insights.

---

## 🌟 Highlights
- **Interactive Dashboard**: The `app.R` script brings all analyses together into a centralized platform, enabling dynamic exploration of results.
- **End-to-End Workflow**: From raw RNA-Seq data to a polished dashboard, this repository covers all stages of data analysis.
- **Custom Gene Insights**: Web scraping ensures detailed descriptions and summaries for all analyzed genes.

---

## 🧑‍💻 Getting Started
1. Clone the repository:
   git clone https://github.com/MarkhamLabVUMC/MarkhamLab.git
2. Follow the installation instructions in the Kallisto_Installation.md and Data_Setup.md files.
3. Explore the dataset by running the notebooks and scripts in the suggested order.

---

## 🚀 Future Enhancements
- 🧬 Add pathway analysis and gene set enrichment functionalities to the dashboard.
- 🎨 Enhance visualizations with additional customization options (e.g., color themes, annotations).
- 🔍 Expand scraping functionality to include pathway or functional data.