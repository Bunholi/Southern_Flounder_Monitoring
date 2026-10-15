# Southern Flounder Monitoring

This README file was updated on 2026-10-10 by Ingrid Vasconcellos Bunholi.

This repository contains the processed data and R scripts used to reproduce the analyses and figures presented in Bunholi et al. *Life-history complexity and population declines constrain molecular monitoring of southern flounder in Texas estuaries*

The analyses combine Texas Parks and Wildlife Department (TPWD) bag-seine, bay-trawl, and gill-net monitoring data with species-specific qPCR and ddPCR results and fish-community DNA metabarcoding results.

## Repository structure

```text
Southern_Flounder_Monitoring_GitHub/
├── README.md
├── data/
│   ├── TPWD_BS_SouthernFlounder_CPUE.csv
│   ├── TPWD_BT_SouthernFlounder_CPUE.csv
│   ├── TPWD_GN_SouthernFlounder_CPUE.csv
│   ├── molecular_occ.csv
│   ├── seine_2024.csv
│   └── trawl_2024.csv
├── scripts/
│   ├── Figure2_TPWD_monitoring.R
│   └── Figures3_4_detection_across_techniques.R
└── outputs/
```

## Data

- `TPWD_BS_SouthernFlounder_CPUE.csv` contains bag-seine catch-per-unit-effort and associated sampling information from 2004–2024.
- `TPWD_BT_SouthernFlounder_CPUE.csv` contains bay-trawl catch-per-unit-effort and associated sampling information from 2004–2024.
- `TPWD_GN_SouthernFlounder_CPUE.csv` contains gill-net catch-per-unit-effort and associated sampling information from 2004–2024.
- `molecular_occ.csv` contains sample-level southern flounder detections from qPCR, ddPCR, and fish-community DNA metabarcoding.
- `seine_2024.csv` and `trawl_2024.csv` contain the 2024 traditional-monitoring observations used for comparison with the molecular methods.
- The LOD/LOQ input dataset contains the standard-curve measurements used with the published Generic qPCR LOD/LOQ Calculator.

## Scripts

### `Figure2_TPWD_monitoring.R`

This script:

- Fits the manuscript's generalized additive models for bag-seine, bay-trawl, and gill-net CPUE trends.
- Fits the binomial detection model used in the manuscript.
- Calculates empirical detection probabilities and confidence intervals for each gear and bay system.

### `Figures3_4_detection_across_techniques.R`

This script:

- Produces the semicircle detection-summary plot comparing traditional fisheries monitoring, metabarcoding, qPCR, and ddPCR (Figure 3).
- Produces the heat map of detection rates by method and bay system (Figure 4).

## Running the analyses

Open R in the repository root, install the required packages if needed, and run:

```r
source("scripts/Figure2_TPWD_monitoring.R")
source("scripts/Figures3_4_detection_across_techniques.R")
```

Required R packages are `dplyr`, `emmeans`, `ggparliament`, `ggplot2`, `mgcv`, `patchwork`, `purrr`, `readr`, `scales`, `stringr`, and `tidyr`.

Generated figures are saved in `outputs/`.

## LOD and LOQ calculation

The limit of detection (LOD) and limit of quantification (LOQ) were calculated using the unmodified **Generic qPCR LOD/LOQ Calculator** developed by Merkes et al. (2019) and described by Klymus et al. (2020). Original calculator: [github.com/cmerkes/qPCR_LOD_Calc](https://github.com/cmerkes/qPCR_LOD_Calc)  

The input dataset used for the LOD/LOQ calculation is included in the `data/` directory. Users can reproduce the calculation by obtaining the original calculator from the links below and running it with that input file.

## References

Klymus, K. E., Merkes, C. M., Allison, M. J., et al. (2020). Reporting the limits of detection and quantification for environmental DNA assays. *Environmental DNA, 2*(3), 271–282. [https://doi.org/10.1002/edn3.29](https://doi.org/10.1002/edn3.29)
Merkes, C. M., Klymus, K. E., Allison, M. J., et al. (2019). *Generic qPCR Limit of Detection (LOD)/Limit of Quantification (LOQ) Calculator* [R script]. U.S. Geological Survey. [https://doi.org/10.5066/P9GT00GB](https://doi.org/10.5066/P9GT00GB)

## Citation

If you use the data or code in this repository, please cite the associated manuscript. The complete manuscript citation will be added here when it becomes available.

## Contact

**Ingrid Vasconcellos Bunholi**  
Email: ingrid.bunholi@utexas.edu  
Affiliation: The University of Texas at Austin Marine Science Institute
