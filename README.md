# Ohio_Bees

Welcome to the Ohio Bees Github folder!



This github involves data analysis of an urban pool derived from compiling data from three papers. The original data from each author is in the folder titled RAW Datasets from authors (see the read me file in said folder).

The project design changed over time, and old iterations of data analysis are stored in the folder titled Scratch. 

The R code file titled Rao function will be needed for both the regional and local analyses.

The file titled Figures contains  1) code used to create some of the figures for the manuscript, 2)a r workspace associated with the figure code, 3) a file titled FiguresforManuscript with the final figures made based on both the regional and local analyses, 4) a file titled Figures_not_used which functions as a scratch folder for the figure creation process. 


Files needed for Regional analysis:

The excel and CSV files labeled bcomm_25.localanalysis contain the presence (1) absence (0) data for the bee species known to the state of Ohio based on the author datasets. 

The excel and csv files labeled btraits_23 list the functional traits associated with each bee species within the regional species pool (state of Ohio). The file titled Categorical Trait key for analysis v3 explains what the functional trait is associated with each number in the btraits_23 files.

The R code file titled Regional to Urban is the code needed for the regional to urban null model analyses comparing the presence of species found as part of the three author data sets to what would be present based on random community assembly. The R workspace titled FinalRegional_workspace is the workspace associated with the Regional to urban code.

The folder titled Final.25.Regional to Urban_Nulls contains the null model csv files based on the regional species pool of Ohio created during the null model creation process.

The folder titled Analysiswithoutseasonalbeestest contains the same r code and regional to urban species pool analysis except with twelve species (Sphecodes aroniae, Sphecodes mandibularis, Andrena neonana, Andrena salictaria, Andrena arabis, Andrena mariae, Pseudopanurgus compositarum, Nomada annulata, Nomada bella, Nomada sulphurate, Epeolus autumnalis, Habropoda laboriosa, Melecta pacifica) removed from the regional species pool. Based on historical collections data these twelve species may not have been active during our sampling period. 



Files needed for local analyses:

The R code titled Final_local_analysis was used to create the csv files titled  traits_urban pool, urbanpool.FS.25, urbanpool.KT.25, and urbanpool.MP.25. These csv files are specific to each author's data set and are sued to run the individual local analyses.

The R code files Final_local_FS, Final_local_KT, and Final_local_MP were used to run the local analyses for each author dataset to compare the bee species found in each habitat type to null models created based on the urban species pool (species found as part of all three studies). Each file has an associated LocalAUTHORINITIALS.workspace file.

The folders titled Urban to Local_Nulls_AUTHORINITIALS contain the null model csv files created during the null model creation process based on the species collected by the three authors. 
