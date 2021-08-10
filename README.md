# aging-omics
# omics pipeline
omics-analysis of gut microbiome and serum metabolome in aging.

# Requirements

- [SOAPnuke](http://manpages.ubuntu.com/manpages/cosmic/man1/soap.1.html)
- [megahit](https://github.com/voutcn/megahit)
- [gmhmmp](http://exon.gatech.edu/license_download.cgi)
- [cd-hit](http://manpages.ubuntu.com/manpages/bionic/man1/cd-hit-para.1.html)
- [NCBI BLAST 2.7.1](https://blast.ncbi.nlm.nih.gov/Blast.cgi)

# A brief description of the contents
Clean and rmhost data/  
quality control for samples including remove low quality reads and remove human sequence  

Assembly/  
Adaptive and iterative assembly for samples with variable coverage.  

Prediction/  
Ab initio gene identification was performed for all assembled contigs.   

Predict gene clustering/  
Creation of non-redundant multi-kingdom skin gene catalog.  

Species annotation/  
Multi-kingdom (bacterial, fungal, viral) taxonomic mapping to NCBI-NT.  

KEGG annotation/  
We aligned putative amino acid sequences, which translated from the ISGC, against the proteins or domains in KEGG databases.  

Quantification of gene, KOand species/  
gene, KO, ARGs and species abundance estimations.  

aging-omics analyses/  
Statistical analysis script and multi-omics analysis script



