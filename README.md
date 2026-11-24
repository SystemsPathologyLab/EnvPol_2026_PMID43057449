# EnvPol_2026_PMID43057449

**[Chronic exposure to low-concentration urban PM2.5 promotes maladaptive kidney repair after ischemic injury via mitochondrial dysfunction and lysosomal stress](https://doi.org/10.1016/j.envpol.2026.129547)**

Sun P, Parra AC, Sanches TR, Wikuats CFH, Butter L, Claessen N, Baelde HJ, Schimmel IM, van der Wel N, Berlingerio SP, Janssens GE, Houtkooper RH, Vaz FM, Roelofs JJTH, Boor P, Strauch M, de Fatima Andrade M, Andrade L, Florquin S, Kers J, Romagnolo A, Tammaro A. *Environ Pollut.* 2026 Nov 19:129547. doi: [10.1016/j.envpol.2026.129547](https://doi.org/10.1016/j.envpol.2026.129547). PMID: [43057449](https://pubmed.ncbi.nlm.nih.gov/43057449/).

For questions about the data or code, please contact Peiqi Sun ([p.sun@amsterdamumc.nl]) or open an issue in this repository. For questions about the paper, please contact the corresponding author, Alessandra Tammaro ([a.tammaro@lumc.nl]).

## Data Availability

Raw RNA-seq data are available in the NCBI Gene Expression Omnibus (GEO) under accession number GSE319410. 
Raw lipidomics and metabolomics data are available in MetaboLights under accession number MTBLS15842.
The processed data matrices and sample metadata for all omics readouts are available in the NCBI Gene Expression Omnibus (GEO) under accession number GSE319410. These can be used as input for the analysis scripts provided in this github. These files allow the downstream analyses and figures to be reproduced without reprocessing the raw experimental data.

## Scripts

### RNAseq.R

This script contains the RNAseq analysis used to generate the following figures:

**Main Figures**

-   Figure 1g – Hallmark pathway enrichment analysis.

-   Figure 3a – Hypoxia-related GO-BP enrichment.

-   Figure 3f – Glycolysis-related GO-BP enrichment.

-   Figure 4b – Lipid catabolic process GSEA.

-   Figure 4f – Differential gene expression volcano plot.

-   Figure 5a – SenMayo GSEA.

-   Figure 7b - Drug repurposing analysis.

**Supplementary Figures**

-   Figure S1b – Hallmark pathway enrichment analysis (Sham+PM vs. Sham+FA).

-   Figure S2a – Hypoxia-related gene expression heatmap.

-   Figure S3a – Fatty acid beta-oxidation GSEA.

-   Figure S3e – GO Cellular Component enrichment analysis.

-   Figure S5a – Response to oxidative stress GSEA.

### Lipidomics_Metabolomics.R

This script contains the lipidomics and metabolomics analysis used to generate the following figures:

**Main Figures**

-   Figure 3e – Metabolomics PLS-DA with top contributing metabolites.

-   Figure 4c – Heatmap of differentially abundant lipid species.

-   Figure 4d – PLSDA of lipid.

-   Figure 4e – Contributors to PLSDA of lipid.

**Supplementary Figures**

-   Figure S2c – Metabolomics volcano plot.

-   Figure S2d – Heatmap of metabolite abundance with metabolic pathway annotations.

-   Figure S3c – Lipidomics volcano plot.

### Pathomics.R

This script contains the pathomics analysis used to generate the following figure:

**Main Figure**

-   Figure 1f – Pathomics PLS-DA.

