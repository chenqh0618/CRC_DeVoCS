# Spatiotemporal Clonal Dynamics of Colorectal Cancer Liver Metastasis Revealed by Dense Volumetric Sampling

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://img.shields.io/badge/DOI-Pending-blue.svg)](https://doi.org/YOUR_PUBLISHED_DOI)

This repository contains the analysis code and key data files to reproduce the findings presented in our manuscript, "Dense Volumetric Crypt-scale Sampling Captures the Spatiotemporal Clonal Dynamics of Colorectal Cancer Liver Metastasis".

## Project Overview

Metastasis is a complex evolutionary process that remains the primary cause of cancer-related mortality. In this study, we employed a dense, three-dimensional, crypt-scale sampling strategy (DeVoCS) to reconstruct the high-resolution evolutionary history of colorectal cancer liver metastasis in two representative cases. Our key findings reveal that:

1.  **Metastasis is an early and continuous process**, initiating long before clinical detection and creating a substantial, clinically occult micrometastatic burden.
2.  **Metastases are founded by multiple, genetically distinct clones (polyclonal seeding)**, suggesting that clonal cooperation may be essential for successful colonization.
3.  **Metastatic tumors evolve along distinct spatial trajectories** following seeding, forming either spatially segregated clonal territories or highly intermixed architectures, a dynamic that appears to be driven by intrinsic cell motility.

This work reframes metastasis not as a singular, late-stage event, but as an early, systemic, and cooperative process, offering new perspectives for therapeutic intervention.


[](figures/analysis_process.png)

---

## Repository Structure

```
.
├── README.md                 # This README file
├── LICENSE                   # MIT License
├── environment.yml           # Conda environment file
│
├── data/                       # Processed input data for analysis scripts
│   ├── metadata/
│   │   └── sample_information.csv
│   ├── genomic_data/
│   │   ├── wgs_somatic_snvs/
│   │   ├── cna_and_purity/
│   │   └── targeted_sequencing/
│   └── reference_data/
│       └── crc_driver_genes.txt
│
├── code/                       # All analysis scripts, organized by workflow
│   ├── 01_wgs_data_processing/
│   ├── 02_heterogeneity_analysis/
│   ├── 03_phylogenetic_reconstruction/
│   ├── 04_clonal_deconvolution_and_seeding/
│   ├── 05_spatial_analysis/
│   ├── 06_molecular_clock_dating/
│   └── utils/
│
├── models/                     # Computational modeling code
│   ├── agent_based_cell_migration/
│   └── continuous_metastatic_seeding/
│
└── results/                    # Generated figures and tables
    ├── figures/
    └── tables/
```

---

## System Requirements & Installation

All analyses were performed on a Linux-based system. We recommend using `conda` to manage the software environment for full reproducibility.

1.  **Clone the repository:**
    ```bash
    git clone https://github.com/chenqh0618/CRC_DeVoCS.git
    cd CRC_DeVoCS
    ```

2.  **Create and activate the conda environment:**
    This will install all required software packages and dependencies (e.g., R, BEAST, IQ-TREE, Python libraries) as specified in the `environment.yml` file.
    ```bash
    conda env create -f environment.yml
    conda activate crc_devocs
    ```

---

## Data Availability

*   **Raw Sequencing Data:** Raw whole-genome and targeted sequencing data have been deposited in the Genome Sequence Archive (GSA) and are publicly available under accession number **PRJCA045427**.
*   **Processed Data:** Key processed data files required to run the analysis scripts in this repository (e.g., filtered VCFs, VAF matrices, sample metadata with spatial coordinates) are available in the `/data` directory.

---

## Reproducibility Workflow

The following steps outline how to reproduce the figures and key results from the manuscript. All scripts should be run from the root directory of the repository (`CRC_DeVoCS/`).

### 1. Heterogeneity Analysis (Results 2.1)

These scripts quantify the multi-scale genetic heterogeneity between and within tumors.

*   **Script:** `code/02_heterogeneity_analysis/2.1_run_pca_and_fst.R`
    *   **Description:** Performs Principal Component Analysis (PCA) and calculates Fixation Index (Fst) on the WGS SNV data.
    *   **Generates:** Figure 1D, Supplementary Figure 3.

*   **Script:** `code/02_heterogeneity_analysis/2.2_run_go_enrichment.R`
    *   **Description:** Conducts Gene Ontology (GO) enrichment analysis on genes with high PCA loadings that separate primary and metastatic samples.
    *   **Generates:** Supplementary Figure 4.

*   **Script:** `code/02_heterogeneity_analysis/2.3_calculate_shannon_diversity.R`
    *   **Description:** Calculates and compares the intra-sample Shannon diversity index using deep targeted sequencing data.
    *   **Generates:** Figures 1E, 1F.

### 2. Phylogenetic and Phylogeographic Reconstruction (Results 2.2)

These scripts reconstruct the evolutionary relationships and trace the geographic migration routes of metastatic seeding.

*   **Script:** `code/03_phylogenetic_reconstruction/3.1_run_iqtree_wgs.sh`
    *   **Description:** Infers a maximum likelihood phylogeny from sparse WGS data using IQ-TREE.
    *   **Generates:** Figure 2A.

*   **Script:** `code/03_phylogenetic_reconstruction/3.2_run_beast_phylogeography.xml`
    *   **Description:** BEAST XML input file for Bayesian phylogenetic and ancestral state reconstruction using the dense targeted sequencing data.
    *   **Generates:** The posterior tree distribution used for Figure 2C.

### 3. Clonal Deconvolution and Seeding Pattern (Results 2.3)

These scripts infer subclonal architectures and test for polyclonal seeding.

*   **Script:** `code/04_clonal_deconvolution_and_seeding/4.1_run_conipher.sh`
    *   **Description:** Runs CONIPHER to deconvolve subclones and reconstruct their phylogeny.
    *   **Generates:** Input for Figures 3A-D.

*   **Script:** `code/04_clonal_deconvolution_and_seeding/4.2_run_treemix.sh`
    *   **Description:** Runs TreeMix to model gene flow between tumor populations, providing independent evidence for polyclonal seeding.
    *   **Generates:** Data for Supplementary Figure 7.

### 4. Spatial Genomics Analysis (Results 2.4)

These scripts analyze the 3D spatial organization of clonal populations.

*   **Script:** `code/05_spatial_analysis/5.1_run_spatial_autocorrelation_and_ibd.R`
    *   **Description:** Performs Moran's I and Isolation-by-Distance (IBD) analyses to quantify spatial structure.
    *   **Generates:** Figures 4A, 4B, Supplementary Figure 8.

*   **Script:** `code/05_spatial_analysis/5.2_generate_3d_clonal_maps.R`
    *   **Description:** Creates 3D visualizations of clonal territories within metastases.
    *   **Generates:** Figures 4D, 4F.

### 5. Molecular Clock Dating (Results 2.5)

These scripts estimate the timing of metastatic dissemination.

*   **Script:** `code/06_molecular_clock_dating/6.2_run_beast_molecular_clock.xml`
    *   **Description:** BEAST XML input file for Bayesian molecular clock analysis to estimate clone divergence times.
    *   **Generates:** The posterior distribution of node ages used in Figures 5A, 5B.

*   **Script:** `code/06_molecular_clock_dating/6.3_plot_divergence_times.R`
    *   **Description:** Processes BEAST output to visualize divergence times and calculate the relative timing of metastasis.
    *   **Generates:** Figures 5A, 5B.

### 6. Computational Modeling

*   **Agent-Based Model of Cell Migration (Results 2.4):**
    *   **Code:** `models/agent_based_cell_migration/`
    *   **Description:** Simulates tumor growth under varying cell migration rates to test its effect on spatial architecture. The script `analyze_abm_results.R` compares simulation output with patient data.
    *   **Generates:** Figure 4G.

*   **Continuous Metastatic Seeding Model (Results 2.5):**
    *   **Code:** `models/continuous_metastatic_seeding/`
    *   **Description:** Implements a mathematical model to simulate the temporal dynamics of metastasis and predict the micrometastatic burden.
    *   **Generates:** Figure 5D, Supplementary Figure 10.

---

## Citation

If you use the code or data from this study, please cite our paper:

> [Your Names], et al. (Year). Title of Your Paper. *Journal Name*, Volume(Issue), pages. [Link to DOI]

**BibTeX:**
```bibtex
@article{YourLastNameYear,
  author    = {[Authors]},
  title     = {Dense Volumetric Crypt-scale Sampling Captures the Spatiotemporal Clonal Dynamics of Colorectal Cancer Liver Metastasis},
  journal   = {[Journal Name]},
  year      = {[Year]},
  volume    = {[Volume]},
  number    = {[Issue]},
  pages     = {[Pages]},
  doi       = {[DOI]}
}
```

## Contact

For questions about the code or analysis, please contact [Your Name] at [your.email@example.com] or open an issue in this repository.

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

