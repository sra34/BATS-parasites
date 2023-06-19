## Parasite-host relationships and links with carbon export in the oligotrophic ocean 

This repository has QIIME 2 compatible files, metadata tables, network files (Cytoscape), and R code needed to perform data analysis for the following manuscript:

**Syndiniales parasites drive species networks and are a biomarker for carbon export in the oligotrophic ocean**<br/>
*Sean R. Anderson, Leocadio Blanco-Bercial, Craig A. Carlson, and Elizabeth L. Harvey,  (2023)*<br/>

* [Preprint](https://www.biorxiv.org)

### 1. Study description 
Protist parasites occupy important roles in marine food webs, influencing host diversity, species succession, and the turnover of organic matter. The marine alveolate group, Syndiniales, are the most diverse and widespread protist parasites in the ocean. Syndiniales have been found in virtually all marine habitats, from surface to deep waters, coastal to offshore, and even in extreme habitats (sediments, oxygen minimum zones, and hydrothermal vents). Yet, their infection dynamics and links with carbon export remain unclear. In this study, we examined a 4-year 18S rRNA metabarcoding dataset from the Bermuda Atlantic Time-series Study (BATS) site in the Sargasso Sea, which included monthly DNA samples collected at twelve depths (1-1000 m). We constructed twelve networks, one for each discrete depth, to observed relationships between Syndiniales and potential hosts. We also leveraged environmental data and POC flux measurements that were collected from the same site. POC flux data were correlated with Syndiniales to provide in situ context for their role in carbon export in marine ecosystems. 

### 2. Bioinformatics processing
* 18S (V4) primers removed with Cutadapt
* Trimmed 18S reads processed with QIIME 2
* Amplicon sequence variants (ASVs) inferred with paired-end DADA2
* Taxonomy assigned for 18S using latest PR2 database release (Version 5.0.1)

### 3. R data analysis and visualization
* QIIME 2 files uploaded to R using [qiime2R](https://github.com/jbisanz/qiime2R)
* Syndiniales population dynamics observed with stacked bar plots, PCoA, and Shannon index
* Environmental drivers of parasites revealed with partial least squares regression (PLSR)
* Covariance networks constructed at each depth using [SPIEC-EASI](https://github.com/zdk123/SpiecEasi)
* Syndiniales-host dynamics observed with depth (number of edges and node degree)
* POC flux at 150 m (sediment traps) correlated to Syndiniales relative abundance 
* Vertical overlap of Syndiniales observed with UpSet plots

### 4. Links to associated data
* Raw sequence data for this project are available in NCBI SRA under BioProject PRJNA769790
* Details on metadata is available at https://bats.bios.edu/data/
