# LOCATOR: anaLysis Of CAncer Tissue micrOenviRonment

We developed a novel computational method, LOCATOR (anaLysis Of CAncer Tissue micrOenviRonment), for spatial analysis of cancer microenvironments using data acquired from mass cytometry imaging (MCI) technologies.  LOCATOR introduces a graph-based representation of tissue images to describe features of the cellular organisation and deploys downstream analysis and visualisation utilities that can be used for data-driven patient risk stratification. Our case study using MCI data from two well-annotated breast cancer cohorts re-confirmed that the spatial organisation of the tumour-immune microenvironment is strongly associated with the clinical outcome in breast cancer. In addition, we report interesting potential associations between the cellular organization of macrophages and patients’ survival.  Our work introduces an automated and versatile analysis framework for MCI data with many applications in future cancer research projects.

# LOCATOR's functions

Before using LOCATOR's functions, some parameters setting are required. Please see the Example.R code in Codes directory for more information. 

Feature matrix can be extracted by:

` > FeatureData <- Extract_Features(Data = MainData)`

The files and plots can be made by main function:

```
> OutPuts <- TIMEClust(FeatureData = FeatureData, 
                       MainData = MainData,
                       SurvivalData = SurvivalData,
                       Cutoff = 'Mean')
```
           
# Example
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

# Contact
Comments and bug reports are welcome, please email: Rezvan Ehsani (rezvan.ehsani@uib.no). We are also interested to know about how you have used our source code, including any improvements that you have implemented. You are free to modify, extend or distribute our source code, as long as our copyright notice remains unchanged and included in its entirety.

# Publication
Article: LOCATOR: feature extraction and spatial analysis of the cancer tissue microenvironment using mass cytometry imaging technologies
DOI: 10.1093/bioadv/vbad146
Journal: Bioinformatics Advances
Pre-print available at bioRxiv
