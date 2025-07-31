# Soil Microbiome Resistance Stability under Environmental Changes

This repository contains R code supporting the research on the impact of community saturation levels on soil microbiome resistance stability, Manuscript in review. 

## Research Overview
**Title**:  
Community saturation levels significantly affect soil microbiome resistance stability under environmental changes  

**Authors**:  
Yuan Chi<sup>a,b,c</sup>, Dong Liu<sup>d</sup>, Zelin Wang<sup>a,b,c</sup>, Kaifang Liu<sup>a,b,c</sup>, Yuan Du<sup>a,b,c</sup>, Yi Duan<sup>a,b,c</sup>, Xin Hu<sup>a,b,c</sup>, Guangbing Xu<sup>a,b,c</sup>, Xinrui Han<sup>a</sup>, Rongxiao Che<sup>a,b,c,*</sup>  
(*Equal contribution: Yuan Chi & Dong Liu*)

**Affiliations**:  
<sup>a</sup> Yunnan Key Laboratory of Soil Erosion Prevention and Green Development, Yunnan University  
<sup>b</sup> State Key Laboratory for Vegetation Structure, Function and Construction, Yunnan University  
<sup>c</sup> Ministry of Education Key Laboratory for Ecosecurity of Southwest China, Yunnan University  
<sup>d</sup> School of Life Sciences, Yunnan University  

**Corresponding Author**:  
Rongxiao Che (cherongxiao@ynu.edu.cn)  

---

### Abstract
Community saturation levels (ratio of observed microbial abundance to environmental carrying capacity) are a critical regulator of microbial community stability. Using a sterilized soil dilution approach across four saturation gradients (10%, 40%, 70%, 100%) in cropland, forest, and grassland ecosystems, we demonstrate that higher saturation levels significantly enhance resistance stability of microbial diversity, community structure, and functional profiles under temperature/moisture/pH/organic carbon/nitrogen manipulations. Community saturation also reshaped co-occurrence networks, increasing robustness while reducing complexity (e.g., 10%-saturation networks had >2× links vs. 100%-saturation). These findings establish community saturation as a key predictor of microbiome stability.

**Keywords**:  
Soil microbiome, Community saturation, Community resistance stability, Environmental changes, Microbial diversity

---

## Data Availability
Raw sequence data are deposited in NCBI SRA under BioProject:  
Raw sequence will be uploaded after manuscript been accepted.
*Amplicon sequencing*: 16S rRNA (bacteria/archaea) & ITS (fungi)  
*Metagenomic data*: Functional profiling  

---

## Code Description
This repository contains **R scripts only** for statistical analyses and visualization.  

### Key Analyses
1. `Amplicon.R`  
   
2. `LDA.R`  
   
3. `Permanova.R`  
   
4. `Metagnomic.R`  


---

## Dependencies
R packages required:  
`tidyverse`, `vegan`, `igraph`, `SpiecEasi`, `lme4`, `phyloseq`, `DESeq2`

---

## Citation
*Manuscript in review* - Please contact authors for citation details.  

---

## Contact
**Corresponding Author**:  
Rongxiao Che  
Yunnan University  
Email: cherongxiao@ynu.edu.cn
