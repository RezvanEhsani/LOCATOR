# LOCATOR: anaLysis Of CAncer Tissue micrOenviRonment


The goal of LOCATOR is to identify subsets of interest cells based on their spatial local microenvironment organization, and find potential clinical and biological associations among the samples which are enriched or decreased for a specific subtype of these interest cells.

Feature matrix can be extracted by:

` > ExtractFeatures.R`

The files and plots can be made by main function:

` OutPuts <- locator(FeatureData = KerenFeatures, 
                   MainData = Data,
                   SurvivalData = SurvivalData,
                   Cutoff = 'Mean')`
                   
An example with downstream analysis can be found in Codes and Example.R file. All returned plots can be found in Keren_Plots directory. The example for Macrophages as cells of interest from Keren et al. paper:

Keren L, Bosse M, Marquez D, Angoshtari R, Jain S, Varma S, Yang SR, Kurian A, Van Valen D, West R, Bendall SC, Angelo M. A Structured Tumor-Immune Microenvironment in Triple Negative Breast Cancer Revealed by Multiplexed Ion Beam Imaging. Cell. 2018 Sep 6;174(6):1373-1387.e19. doi: 10.1016/j.cell.2018.08.039. PMID: 30193111; PMCID: PMC6132072.
                   
                   
