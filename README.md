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

The information about the original (mzML) files is available in the files:

- [get_rt.Rmd](get_rt.Rmd): Markdown file for obtaining the RTs of the 
  maximum peak corresponding to each compound (considering its exact mass). 
  This serve as a starting point for define the RT of each compound.  
  - `study`: specify which of the 2 studies we want to focus 
  (i.e., "internal_standards" or "standards_dilution").
  - `mixnum`: in case we specified `study <- "standards_dilution"`, 
  specify to which MIX we want to focus.
  - `polarity`: specify in which polarity ("POS" or "NEG") we want to focus.
- [EIC_superposed.Rmd](EIC_superposed.Rmd): 
- [EIC_manually.R](EIC_manually.R): R file for plot the EIC for a 
  specific compound in both ionization modes.
- [spectras_max_ions.Rmd](spectras_max_ions.Rmd): 
- [spectra_manually.R](spectra_manually.R): 
- [internal_standards_evaluation.Rmd](internal_standards_evaluation.Rmd): 
- [standard_dilution_dose_response.Rmd](standard_dilution_dose_response.Rmd): 
- [exclusion_mz.R](exclusion_mz.R): 
