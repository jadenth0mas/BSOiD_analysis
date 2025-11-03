# Seizure Behavior Decoding and SUDEP Risk Prediction

**Statistical analysis code for quantifying seizure microfeatures and predicting sudden death risk in epilepsy mouse models**

[![Publication](https://img.shields.io/badge/Publication-Annals%20of%20Neurology-blue)](https://doi.org/10.1002/ana.78032)

## Overview

This repository contains the statistical analysis pipeline for behavior decoding in mouse models of epilepsy. The project applies multivariate statistical methods to high-dimensional behavioral data to:

- Delineate seizure states across different genotypes
- Identify microfeatures associated with Sudden Unexpected Death in Epilepsy (SUDEP)
- Model seizure progression over time using longitudinal data

This work was conducted at The Ohio State University Gu Lab and contributed to a publication in *Annals of Neurology*.

## Research Context

**Problem:** SUDEP is a leading cause of epilepsy-related mortality, but predicting risk remains challenging.

**Approach:** We used quantitative behavioral analysis combined with advanced statistical modeling to identify seizure characteristics predictive of sudden death risk.

**Impact:** Our methods enable objective, data-driven seizure classification and risk stratification.

## Methods

### Full Analysis Pipeline

This project integrates computer vision, unsupervised learning, and classical statistics:

```
Video Data → DeepLabCut → B-SOID → Statistical Analysis → SUDEP Risk Models
(pose estimation) (behavior clustering) (this repository)
```

#### Upstream Processing (Not in this repo)
- **DeepLabCut**: Markerless pose estimation to extract animal keypoint coordinates from video
- **B-SOID (Behavioral Segmentation of Open-field in DeepLabCut)**: Unsupervised learning to identify behavioral motifs and generate time-series expressions of behavior states

#### Statistical Analysis (This repo)
- **Principal Component Analysis (PCA)**: Dimensionality reduction of behavioral features and B-SOID expression patterns
- **Linear Discriminant Analysis (LDA)**: Seizure state classification across genotypes using behavior group features
- **Hierarchical Clustering**: Identification of behavioral phenotype groups from B-SOID outputs
- **Mixed Effects Models**: Longitudinal analysis of seizure progression and behavior expression over time

### Analysis Pipeline
1. **Data preprocessing**: Import B-SOID behavioral expression data and seizure annotations
2. **Feature engineering**: Extract summary statistics from behavior group time-series (duration, frequency, transitions)
3. **Exploratory analysis**: PCA for visualization of behavior space and pattern identification
4. **Classification**: LDA to distinguish seizure states and genotype-specific patterns using behavior features
5. **Clustering**: Hierarchical methods to identify natural groupings based on behavioral expression profiles
6. **Longitudinal modeling**: Mixed effects models for temporal progression of behavior patterns and SUDEP risk

## Key Findings

- Successfully delineated seizure microfeatures with high discriminative power
- Identified behavioral patterns associated with increased SUDEP risk
- Demonstrated genotype-specific seizure progression patterns
- Validated classification models across multiple mouse strains

## Technologies

**Language:** R 4.4.1

**Key Packages:**
- `tidyverse`: Data manipulation and visualization
- `lme4`: Mixed effects modeling
- `MASS`: LDA implementation
- : PCA visualization
- `cluster`: Hierarchical clustering
- `ggplot2`: Publication-quality figures

**Related Tools (upstream processing, not included):**
- [DeepLabCut](https://github.com/DeepLabCut/DeepLabCut): Pose estimation from video
- [B-SOID](https://github.com/YttriLab/B-SOID): Unsupervised behavioral segmentation

## Data Availability

**Note:** The behavioral data used in this analysis is not uploaded in this repository. The code is provided to demonstrate statistical methodology and analytical approach. 

Researchers interested in the data should contact the Gu Lab at The Ohio State University.

## Installation & Usage

### Prerequisites
```r
# Install required packages
install.packages(c("tidyverse", "lme4", "MASS", "factoextra", 
                   "cluster", "ggplot2", "cowplot"))
```

### Running the Analysis

The scripts are numbered in the order they should be run. Each script is self-contained with detailed comments.

```r
# Example: Run PCA analysis
source("02_exploratory_pca.R")

# Example: Fit mixed effects model
source("05_mixed_effects_models.R")
```

### Expected Inputs

Each analysis script expects data derived from the B-SOID pipeline:

**Behavioral Expression Data:**
- Time-series of behavior group assignments (output from B-SOID)
- Behavioral bout features: duration, frequency, transition probabilities
- Pose-derived metrics from DeepLabCut keypoints (when applicable)

**Metadata:**
- Subject ID, genotype, experimental timepoint
- Seizure annotations and outcomes (including SUDEP events)
- Experimental conditions and covariates

**Data Structure Example:**
```r
# Expected format for behavioral expression
behavior_data <- data.frame(
  subject_id = c("M001", "M001", ...),
  timepoint = c(1, 2, ...),
  behavior_group = c(3, 7, 3, ...),  # B-SOID cluster assignment
  bout_duration = c(2.3, 1.5, ...),
  genotype = c("WT", "WT", ...),
  sudep_outcome = c(0, 0, ...)
)
```

See inline comments in each script for specific requirements.

## Publication

This work is published in *Annals of Neurology*:

**Shen, Y., Thomas, J., Chen, X., Zelidon, J., Najeeb, M., Hahn, A., Zhang, P., Sathyanesan, A. and Gu, B. (2025).** Behavior Decoding Delineates Seizure Microfeatures and Associated Sudden Death Risks in Mouse Models of Epilepsy. *Annals of Neurology*. [https://doi.org/10.1002/ana.78032](https://doi.org/10.1002/ana.78032)

## Author

**Jaden Thomas**  
M.S. Biostatistics Student, University of Minnesota  
B.S. Data Analytics, The Ohio State University

Statistical analysis conducted as Student Research Assistant in the Gu Lab (March 2024 - September 2025)

## Acknowledgments

- **Principal Investigator:** Dr. Bin Gu, The Ohio State University
- **Collaborators:** Gu Lab members and co-authors
- Funded by [include if you know the funding source]

## License

Code is provided for educational and research purposes. Please cite the publication if you use or adapt these methods.

---

*For questions about the methodology or code, please open an issue or contact me via [your email/LinkedIn].*
