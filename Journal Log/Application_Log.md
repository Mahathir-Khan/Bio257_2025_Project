# 📓 Application Log

This log tracks updates and changes made to the repository. Each entry corresponds to one or more commits.

---
TEMPLATE: DO NOT DELETE

## 🗓️ [Date: 2025-10-03]
**Name:** [Mahathir] 
**Description:**  
- Github repo was established
- Folders for Papers, Data, Application Log, and ReadMe established
  
---

## 🗓️ [Date: 2025-10-07]
**Name:** [Mahathir] 
**Description:**  
- Added Group Agreement pdf to front page of repo
- Added Project outline to repo main thread
**Name:** [Aman] 
**Description:**  
- Created template of outline
- Hunted for literature and relevant datasets


---
## 🗓️ [Date: 2025-10-14]
**Name:** [Sofia Kaplan] 
**Description:**  
- [I added data files for the 16s (Emery et al. 2017) and single cell RNA seq data and determined which studies we will be using (Darmanis et al. 2015)]

---

## 🗓️ [Date: 2025-10-16]

**Name:** [Aman] 
**Description:**  
- Read papers (2) of dataset
    - Report findings
- Dig through for pac, umap, tsne
**Name:** [Kelechi] 
**Description:**  
- Made goals for our meeting on sunday
- Litertaure Review.

## 🗓️ [Date: 2025-10-18]
**Name:** [Sofia] 
**Description:**  
- Updated single cell RNA seq paper to Wang et al 2024 and added databases to Papers for Data (corrected to include AD vs control data)

**Name:** [Mahathir]
**Description:**
- Set up the project notebook and am reading up on literature for the phylogenetic analysis we
intend to perform.
- Hope to have specific libraries and methods ready by next week

## 🗓️ [Date: 2025-10-19]
**Name:** [Mahahthir] 
**Description:**  
- Found papers for PhyloSeq Workflow. Now identifying interesting graphics to generate and intends to complete phylogenetic analysis within 2 weeks.
- graphics of interest are: histogram of dominant taxa, heatmap of samples vs taxa, phylogenetic tree

**Name:** [Mahahthir, Aman, Sofia, Kaycee] 
**Description:**  
- Had weekly team meeting to discuss project updates and perform a literature review
  
---

## 🗓️ [Date: 2025-10-24]
**Name:** [Sofia] 
**Description:**  
- Found new paper for RNA-seq data Maffioli E et al. 2022.
- Downloaded data to bluehive and converted as fastq files.

**Name:** [Mahathir]
**Description:**
- contributed to the outline
- Helped search for new datasets
- Accumulated accession numbers for RNAseq data to be uploaded to bluehive
- Read into statistical tests
- Swapped analysis priorities with Aman. Moving forward, I will focus on bulk RNAseq analysis with Sofia, while Aman and Kalechi focus on 16s phylogenetic analysis.

## 🗓️ [Date: 2025-10-23]
**Name:** [Aman] 
**Description:**  
- Began looking for pipelines for 16s RNA and bulk RNA-seq data
- Outlined Presentation 
---

## 🗓️ [Date: 2025-10-31]
**Name:** [Aman] 
**Description:**  
- Constructed DADA2 pipeline as an alternative to 16s RNA fastq file pre-processing

**Name:** [Mahathir]
**Description:**
- Finalized reference genome selection
- Begin writing trimming script for paired-end reads. Will consult with TAs regarding appropriate memory usage before running the script.
- Presented proposal to class and did my portion of the slide deck

---
## 🗓️ [Date: 2025-11-7]
**Name:** [Aman] 
**Description:**  
- Redownloaded dataset and adjusting pipeline due to initial formatting mishaps with 16s data
  - Downloaded incorrect reverse fastq sequences as database did not clearly indicate forward vs reverse reads
  - Database also updated so only contains 16s data
    - Concerned may have been error with previous data, so redownloaded just in case
- Downloaded bash script from database with series of wget commands to import fastq files
  - wget commands incorrectly written so corrected locally before sending file to bluehive
  - ran bash script to import all files to bluehive

**Name:** [Mahathir]
**Description:**
- Clean github to consolidate all datasets and accession numbers into one folder in data directory titled Final_Data_Selections
- Clean through references and compile in papers fold in md file titled References_master_list
- Complete script for trimming RNAseq data at phred 20
- Script took roughly 6 hours to run using 32 gb memory and outputted trimmed files into Bluehive directory labeled trimmed_fastq_for_RNAseq
- Awaiting quality check of trimmed RNA by Sofia and will then proceed with
pseudomapping to transcriptome with Salmon
- Cleaned up Application Log and submitted Journal Discussion Questions 
---
## 🗓️ [Date: 2025-11-9]
**Name:** [Sofia] 
**Description:**  
 - Finished generating FastQC reports and created a table of the trimming reports. File uploaded under Data as https://github.com/Mahathir-Khan/Bio257_2025_Project/blob/main/Data/trimming_summary_phred20.tsv
   
**Name:** [Aman] 
**Description:**  
- Corrected pipeline to trim and filter fastq data
  - Removed Chimeras
  - Calculated error rates across reads
- Processed files into ASV format to prepare for additional processing with vegan and phyloseq
---

## 🗓️ [Date: 2025-11-14]
**Name:** [Sofia] 
**Description:**  
 - Uploaded ssREAD data to bluehive under /scratch/bio257_2025/Users/group4_ADxTummy/RNAseq_data/ssread/ssread_files.
---

**Name:** [Aman] 
**Description:**  
- Helped Kaycee to prepare 16s data analysis pipeline
- ASV did not contain arm information (control or AD)
  - Manually created tsv file containing arm assignments associated with each accession number
  - Merged arm assignment info into ASV
---

**Name:** [Mahathir] 
**Description:**  
- Downloaded human reference assembly into bluehive, and converted to fa.gz file
- Organized bluehive file structure, creating new directories for salmon outputs
- Begin writing script for Salmon transcriptome alignment

## 🗓️ [Date: 2025-11-18]
**Name:** [Aman] 
**Description:**  
- Downloaded git directory associated with ssread data processing
---

## 🗓️ [Date: 2025-11-23]
**Name:** [Mahathir] 
**Description:**  
- Generated salmon index using salmon index.sh (in code directory)
- Ran salmon script and obtain quant files for all 12 samples with salmon_script.sh (in code directory)
- Proceeded with differential expression analysis of quant information in Salmon_Analysis.Rmd (in code directory)
- Produced downregulated_genes_AD_vs_Control.csv (in data folder)
- Produced upregulated_genes_AD_vs_Control.csv (in data folder)
---

- ## 🗓️ [Date: 2025-11-18]
**Name:** [Aman] 
**Description:**  
- Assisted Sofia in debugging ssRead pipeline
- Began outlining poster
