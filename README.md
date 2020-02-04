# Measurement of standards in different LC/MS configurations

This repository contains measurement of sets of standards on LC-MS systems with
specific configurations. The aim is to provide approximate retention times 
(RTs) for the standard compounds in the respective setups.  
  

## Data tables

- [internal_standards_files.txt](data/internal_standards_files.txt): table 
  indicating which mzML files correspond to the experiment in which 
  internal standards (solved in water) were repeatedly measured 
  (from September 2019).
- [internal_standards.txt](data/internal_standards.txt): table containing 
  formulas (or neutral exact masses), highest ions in each mode (POS & NEG), 
  and retention times for each standard.
- [standards_dilution_files.txt](data/standards_dilution_files.txt): table 
  indicating which mzML files correspond to the experiment `Dynamic Range` 
  in which sets of standards (solved in water) were measured (in different
  concentrations).
- [standards_dilution.txt](data/standards_dilution.txt): table containing 
  formulas (or neutral exact masses), highest ions in each mode (POS & NEG), 
  and retention times for each standard.
- [exclusion_mz.txt](data/exclusion_mz.txt): table containing the 
  (exact) masses of mz values that probaly are noise / backgroun ions. 
  This is a "dynamic" table that is being updated along the time.


## R files

- [get_rt.Rmd](get_rt.Rmd): Markdown file for obtaining the RTs of the 
  maximum peak corresponding to each compound (considering its exact mass). 
  This serve as a starting point for define the RT of each compound. 
  At the begining we have to specify the following parameters:  
  - `study`: specify which of the 2 studies we want to focus 
  (i.e., "internal_standards" or "standards_dilution").
  - `mixnum`: in case we specified `study <- "standards_dilution"`, 
  specify to which MIX we want to focus.
  - `polarity`: specify in which polarity ("POS" or "NEG") we want to focus.
- [EIC_superposed.Rmd](EIC_superposed.Rmd): Markdown file for printing 
  the EIC of one compound using all samples. The plots are saved 
  in the folder `plots`.
- [EIC_manually.R](EIC_manually.R): R file for plot the EIC for a 
  specific compound in both ionization modes.
- [spectras_max_ions.Rmd](spectras_max_ions.Rmd): Markdown for printing 
  the "cleaned" spectras of all compounds in a sample (i.e., after excluding 
  the mz values specified in the file `exclusion_mz.txt`)
- [spectra_manually.R](spectra_manually.R): R file for plot the spectra 
  for a specific compound.
- [internal_standards_evaluation.Rmd](internal_standards_evaluation.Rmd): 
  Markdown file that plot a graph with 1 boxplot / injection and all its  
  measures, and that generates a pdf file with the trends of each compound 
  along the samples. 
- [standard_dilution_dose_response.Rmd](standard_dilution_dose_response.Rmd): 
  Markdown file that generates a pdf file showing the linearity of 
  each compound.
- [exclusion_mz.R](exclusion_mz.R): R code for generate the 
  table `exclusion_mz.txt`.

## Matrix effects
In the folder `std_serum` there are the files used for analyse the 
experiments about IS and STDs (injected in different MIXs) 
at 2 different concentrations in water and in QC-serum samples.  
The code [RT_matrix_effect](std_serum/RT_matrix_effect.Rmd) plots all 
the EICs present in the sample (ie, IS or MIX 01-20) in the folder `images`.  
In that folder there is also the code [`MS2_plot`](std_serum/MS2_plot.Rmd) 
to search a MS2 spectrum for a specific mz-rt value

