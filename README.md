# Dormitory Culture Analysis

[![DOI](https://img.shields.io/badge/DOI-10.31235/osf.io/t3nbf_v1%2FXXXXX-blue)]([https://doi.org/10.31235/osf.io/t3nbf_v1](https://doi.org/10.31235/osf.io/t3nbf_v1))
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

## Overview

This repository contains the analysis code and manuscript for the study:

**"Residential Configuration and Dormitory Culture Formation: A Quantitative Analysis of How Gender Composition and Spatial Design Influence Residents' Cooperative Attitudes in University Housing"**

Authors: Tsubasa Sato, Natsumi Negoro, Misa LoPresti, Misaki Iio  
Affiliation: Keio University

## Abstract

This study examined how different gender compositions (male-only, female-only, and co-ed dormitories) and spatial designs (same-floor co-ed and floor-separated co-ed) across four buildings in the H village student dormitory at Keio University influenced residents' cooperative attitudes and values. A longitudinal survey was conducted five times between May and November 2025, and data from 153 domestic students were analyzed.

Key findings:
- The overall ICC was 0.075 (95% CI: [0.000, 0.258]), indicating that approximately 7.5% of response variance was explained by building differences
- For "willingness to participate in dormitory events," ICC reached 0.279 (p < 0.001)
- No significant difference was found between similarity to same-sex seniors (95.1%) and opposite-sex seniors (95.7%; p = 0.160)
- Co-ed buildings showed significantly lower event participation willingness compared to single-sex buildings (p < 0.05)

## Repository Structure

```
dormitory-culture-analysis/
├── README.md
├── LICENSE
├── CITATION.cff
├── .gitignore
├── R/
│   ├── 01_data_preprocessing.R
│   ├── 02_icc_analysis.R
│   ├── 03_cultural_similarity.R
│   ├── 04_gender_influence.R
│   ├── 05_building_comparison.R
│   ├── 06_visualization.R
│   └── analysis_main.R
├── manuscript/
│   ├── h_village_paper_socarxiv.tex
│   ├── h_village_paper_japanese.tex
│   └── references.bib
├── figures/
│   ├── fig2_cultural_effects.png
│   ├── fig3_convergence.png
│   ├── fig4_gender_influence_detail.png
│   ├── fig5_building_type_comparison.png
│   ├── fig6_q2_barrier.png
│   └── fig7_building_type_overview.png
└── data/
    └── README.md (data access information)
```

## Requirements

### R Packages

```r
install.packages(c(
  "lme4",
  "vegan",
  "boot",
  "ggplot2",
  "dplyr",
  "tidyr",
  "readr",
  "patchwork",
  "ggpubr",
  "RColorBrewer"
))
```

### Software Versions

- R version 4.3.0 or later
- RStudio (recommended)

## Usage

### Running the Analysis

1. Clone this repository:
```bash
git clone https://github.com/Rhizobium-gits/dormitoly-culture-analysis.git
cd dormitoly-culture-analysis
```

2. Open RStudio and set the working directory to the repository root.

3. Run the main analysis script:
```r
source("R/analysis_main.R")
```

### Compiling the Manuscript

The manuscript is written in LaTeX for SocArXiv submission:

```bash
cd manuscript
pdflatex h_village_paper_socarxiv.tex
bibtex h_village_paper_socarxiv
pdflatex h_village_paper_socarxiv.tex
pdflatex h_village_paper_socarxiv.tex
```

Note: The SocArXiv template (`sa.cls`) is required for compilation.

## Data Availability

Due to privacy considerations for survey participants, the raw data cannot be made publicly available. Requests for data access for research purposes should be directed to the corresponding author at satotsubasa@keio.jp.

## Statistical Methods

### Intraclass Correlation Coefficient (ICC)

Building cultural influence was estimated using linear mixed-effects models:

$$y_{ij} = \beta_0 + u_j + \epsilon_{ij}$$

where:
- $y_{ij}$ is the score of individual $i$ in building $j$
- $\beta_0$ is the overall mean
- $u_j \sim N(0, \sigma^2_b)$ is the building random effect
- $\epsilon_{ij} \sim N(0, \sigma^2_e)$ is the residual

ICC was calculated as:

$$\text{ICC} = \frac{\sigma^2_b}{\sigma^2_b + \sigma^2_e}$$

### Cultural Similarity

Cosine similarity between individual and group profiles:

$$\text{Similarity} = \frac{\sum_{i=1}^{8} x_i y_i}{\sqrt{\sum_{i=1}^{8} x_i^2} \sqrt{\sum_{i=1}^{8} y_i^2}} \times 100\%$$

### PERMANOVA

Permutational multivariate analysis of variance based on Aitchison distance was used to test multivariate differences between buildings.

## Citation

If you use this code or reference this study, please cite:

```bibtex
@article{sato2025dormitory,
  title={Residential Configuration and Dormitory Culture Formation: A Quantitative Analysis of How Gender Composition and Spatial Design Influence Residents' Cooperative Attitudes in University Housing},
  author={Sato, Tsubasa and Negoro, Natsumi and LoPresti, Misa and Iio, Misaki},
  journal={SocArXiv Preprint},
  year={2025},
  doi={10.31235/osf.io/XXXXX}
}
```

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## Contact

- Tsubasa Sato - okay.bio.sato at gmail.com
- Misaki Iio - misaki.iio at keio.jp

## Acknowledgements

We thank the H village Resident Assistants (RAs) and all residents who participated in this survey. We also thank the H village management office for their cooperation in conducting this study.
