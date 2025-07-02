# swiss_co-circulating_viruses
Associated code with for the study of respiratory viruses sequenced in Switzerland during the 2023/24 season.

The repository is structured as follows:
```
.
├── Data
├── Overview_circulating_viruses
│   ├── fig_1_data_format.R
│   └── fig_1_plot.R
├── Phylogenetic_analysis
│   └── README.md
├── README.md
├── Sequencing_results
│   ├── README.md
│   ├── fig_2_data_format.R
│   ├── fig_2_plot.R
│   ├── fig_3_data_format.R
│   └── fig_3_plot.R
├── Supplements
└── analysis_pipeline.Rmd
```

- Data repository contains all pertinent data for the analysis
- Overview_circulating_viruses directory contains scripts to generate the first figure
- Phylogenetic_analysis directory contains the Nextstrain pipelines to generate all phylogenies and the scripts to generate the tree figures of the manuscript
- Sequencing_results directory contains the scripts to generate the figures pertaining to the sequencing results, including coverage plots
- Supplements contains the scripts to generate the supplementary figures and the coverage plots of all single-stranded viruses from this study
- analysis_pipeline.Rmd contains the main analysis of this study 

# Running the pipeline

You must unzip the depth files (Data/depth_files/) to obtain the coverage figures. 
