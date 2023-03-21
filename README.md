# LOCATOR: anaLysis Of CAncer Tissue micrOenviRonment


The goal of LOCATOR is to identify subsets of interest cells based on their spatial local microenvironment organization, and find potential clinical and biological associations among the samples which are enriched or decreased for a specific subtype of these interest cells.

Feature matrix can be extracted by:

` > ExtractFeatures.R`

The files and plots can be made by main function:

` 
> OutPuts <- TIMEClust(FeatureData = KerenFeatures, 

                   MainData = Data,
                   
                   SurvivalData = SurvivalData,
                   
                   Cutoff = 'Mean')
                   `
                   
An example can be found in Codes directory and Example.R file. All returned plots can be found in Keren_Plots directory. The example is for Macrophages as cells of interest from Keren et al. paper:

Keren L, Bosse M, Marquez D, Angoshtari R, Jain S, Varma S, Yang SR, Kurian A, Van Valen D, West R, Bendall SC, Angelo M. A Structured Tumor-Immune Microenvironment in Triple Negative Breast Cancer Revealed by Multiplexed Ion Beam Imaging. Cell. 2018 Sep 6;174(6):1373-1387.e19. doi: 10.1016/j.cell.2018.08.039. PMID: 30193111; PMCID: PMC6132072.
                   
# Dependencies and System Requirements
LOCATOR is developed with R 4.2.2 (2022-10-31 ucrt).

The analysis presented in our publication requires packages to be installed, and we refer to the following document about the R session info.

` mySessionConfig.sh`

# License
This project is licensed under the MIT License.

Copyright 2022 Department of Informatics, University of Bergen (UiB) and the Centre of Cancer Biomarkers (CCBIO), Norway

You may only use the source code in this repository in compliance with the license provided in this repository. For more details, please refer to the file named "LICENSE.md".
