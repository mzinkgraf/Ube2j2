# Ube2j2
This package implements the analysis of mass spectrometry data from Ube2j2 deficient cell lines using a WGCNA network approach.

Ube2j2 can be installed in R using the following commands:

````{r}
    #install dependencies
    install.packages(c("WGCNA", "reshape2", "fastcluster", "dynamicTreeCut", 
        "flashClust", "plyr", "ggplot2", "gridExtra", "ggdendro", "emmeans", "stringr"), 
        dependencies = TRUE)

    #install Ube2j2
    require(devtools);
    install_github("mzinkgraf/Ube2j2");
    
````

The R script containing the analysis can be found in `R/Ube2j2_analysis.R` and supporting files are contained in the `Data/` folder.

Plots of individual genes and how they respond to experimental treatments can be found in `Data/plots/`

> [!IMPORTANT]
> Note on plots: Github only shows the first 1000 proteins in the folder and others can be obtained using the `Go to file` search box in upper left with the uniprot name.
<br/><br/>
<img src="Data/Go2file.png"  width="50%">
