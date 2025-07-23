# swiss co-circulating viruses

Associated code with for the study of respiratory viruses sequenced in Switzerland during the 2023/24 season, manuscript called: "Characterizing Co-Circulating Respiratory Virus Genomic Diversity in Switzerland with Hybrid-Capture Sequencing and Phylogenetic Reconstructions: Insights into the 2023/24 Season". 

## Repository structure

The repository is structured as follows:

```
.
├── analysis_pipeline.md
├── analysis_pipeline.Rmd
├── coverage_plots_all.md
├── coverage_plots_all.Rmd
├── Data
│   ├── data
│   └── resources
├── Figures
├── Overview_circulating_viruses
│   ├── fig_1_data_format.R
│   └── fig_1_plot.R
├── Phylogenetic_analysis
│   ├── Nextstrain_builds
│   ├── README.md
│   └── tree_plot.R
├── README.md
├── Sequencing_results
│   ├── coverage_genome_data_format.R
│   ├── coverage_genome_plot.R
│   ├── fig_2_data_format.R
│   ├── fig_2_plot.R
│   ├── fig_3_data_format.R
│   ├── fig_3_plot.R
│   └── README.md
└── Supplements
    ├── fig_s2_data_format.R
    ├── fig_s2_plot.R
    └── README.md
```

markdown documents include: `analysis_pipeline.md` and `coverage_plots_all.md`. `analysis_pipeline.md` contains all the compiled main figures and code used in the main manuscript and `coverage_plots_all.md` contains the coverage plots of all single-stranded viruses from this study. 

* `analysis_pipeline.Rmd` contains the code of the main analysis of this study
* `coverage_plots_all.Rmd` contains the code for the coverage plots of all single-stranded viruses from this study
* `Data/` contains all pertinent data for the analysis and the resources needed to run the pipeline
* `Overview_circulating_viruses/` contains scripts to generate figure 1
* `Phylogenetic_analysis/` contains (1) the Nextstrain pipelines to generate all phylogenies and (2) the scripts to generate the tree figures of the manuscript
* `Sequencing_results/` contains the scripts to generate the figures pertaining to the sequencing results, including coverage plots (figures 2,3,4)
* `Supplements/` contains the scripts to generate the supplementary figures and tables.

## Running the pipeline

You can reproduce the results from this study by running `analysis_pipeline.Rmd`, with required packages installed and detailed in the file. However, you must populate the depth file folder (Data/depth_files/) to run the sequencing depth analysis. This data is available upon request.

To run Nextstrain builds, head to `Phylogenetic_analysis/Nextstrain_builds/` and into individual virus subdirectories to run every pipleine. A detailed ReadMe file is incldued in each subdirectory. 

Interactive Nextstrain visualizations are available [here](https://nextstrain.org/community/cevo-public/ReVSeq-project)