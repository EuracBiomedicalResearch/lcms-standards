# Measurement of standards in different LC/MS configurations

**version 1.0.0.0** (2025-01-15)

This repository contains measurement of sets of standards on LC-MS systems with
specific configurations. The aim is to provide approximate retention times 
(RTs) and ions for the standard compounds in the respective setups.

The first definition of retention times and ions for all standards was performed
manually by Mar Garcia-Alloy. In the *new* workflow we aim to use the `xcms` and
`Spectra` package to refine these and to determine also all additional adducts
created/measurable from each compound. In the workflow we first:

- Identify chromatographic peaks in each file for one *standard mix*.
- Perform correspondence analysis to group each chromatographic peak to a
  feature.
- Identify features with matching *m/z* and approximate retention time.
- Plot and (visually) inspect their EICs.
- Extract their MS2 spectra and match them against HMDB and MassBank *reference*
  spectra.
- For validated features, add their information to the *ion database*.

## Annotation database

The annotation results of the standards is provided in the SQLite database
*IonDb.IfB_HILIC.1.0.0.0.sqlite*. This database is in the *IonDb* format of the
[CompoundDb](https://github.com/RforMassSpectrometry/CompoundDb) Bioconductor
package and it can be easily integrated into annotation workflows using the
[MetaboAnnotation](https://github.com/RforMassSpectrometry/MetaboAnnotation)
Bioconductor package. The database provides retention times and *m/z* values for
ions of the analyzed standards in database table columns *ion_rt* and
*ion_exactmz*.

Below is an example how the database could be loaded and how relevant
information for annotation with *MetaboAnnotation* could be extracted.

```
library(CompoundDb)
library(MetaboAnnotation)

#' Load the data as an IonDb database
idb <- IonDb("IonDb.IfB_HILIC.1.0.0.0.sqlite")

#' Extract ion information from the database
std_ions <- ions(idb, 
                 columns = c("compound_id", "name", "ion_adduct", 
				             "ion_rt", "ion_exactmz", "polarity"))

#' To annotate compounds for positive polarity
std_ions <- std_ions[std_ions$polarity == "POS", ]

#' Assuming lcms_features is a data.frame with data for LC-MS features
#' from an LC-MS analysis, e.g. using the xcms package. The columns
#' mzmed and rtmed in that table would contain the m/z and retention
#' time values for each feature. Annotation could then be done for
#' example with
ann <- matchValues(
    lcms_features, std_ions, MzRtParam(ppm = 20, toleranceRt = 4),
	rtColname = c("rtmed", "ion_rt"), mzColname = c("mzmed", "ion_exactmz"))
```


## Analysis workflow files

- [match-standards-introduction.Rmd](match-standards-introduction.Rmd): A 
  detailed description of the approach per *standard mix* is given.
- [match-standards-serum-mix01.Rmd](match-standards-mix01.Rmd): matching and
  identifying signal from standards of mix 01. This includes preprocessing and
  ultimately defines the retention time, the measured ions and related MS/MS
  spectra of the standards in serum.
- ...
- [match-standards-serum-mix20.Rmd](match-standards-mix08.Rmd): matching and
  identifying signal from standards of mix 08. This includes preprocessing and
  ultimately defines the retention time, the measured ions and related MS/MS
  spectra of the standards in serum.
- [create_reference_database.Rmd](create_reference_database.Rmd): document
  describing the generation of the reference database that can be used for
  annotation of LC-MS experiments that used the IfB HILIC HPLC-MS setup.
  

### Standard mixes process and status

#### Standards spiked with QC Pool samples

- [X] mix 01. 2024-12. Created by Andrea, checked by Jo.
- [X] mix 02. 2024-12. Created by Andrea, checked by Jo.
- [X] mix 03. 2025-01. Created by Andrea, checked by Jo.
- [X] mix 04. 2025-01. Created by Andrea, checked by Jo.
- [X] mix 05. 2025-01. Created by Andrea, checked by Jo.
- [X] mix 06. 2025-01. Created by Andrea, checked by Jo.
- [X] mix 07. 2025-01. Created by Andrea, checked by Jo.
- [X] mix 08. 2025-01. Created by Andrea, checked by Jo.
- [X] mix 09. 2024-11. Created by Andrea, checked by Jo.
- [X] mix 10. 2024-11. Created by Andrea, checked by Jo.
- [X] mix 11. 2024-11. Created by Andrea, checked by Jo.
- [X] mix 12. 2024-12. Created by Marilyn, checked by Jo.
- [X] mix 13. 2024-12. Created by Marilyn, checked by Jo.
- [X] mix 14. 2024-12. Created by Marilyn, checked by Jo.
- [X] mix 15. 2024-12. Created by Marilyn, checked by Jo.
- [X] mix 16. 2024-12. Created by Marilyn, checked by Jo.
- [X] mix 17. 2024-12. Created by Marilyn, checked by Jo.
- [X] mix 18. 2024-12. Created by Marilyn, checked by Jo.
- [X] mix 19. 2024-12. Created by Marilyn, checked by Jo.
- [X] mix 20. 2024-12. Created by Marilyn, checked by Jo.


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
  and retention times for each standard. This file is now obsolete. It was used
  as basis for the new annotation and adduct identification described above. For
  annotation, the new *IonDb.IfB_HILIC* SQLite database should be used.
- [exclusion_mz.txt](data/exclusion_mz.txt): table containing the 
  (exact) masses of mz values that probaly are noise / backgroun ions. 
  This is a "dynamic" table that is being updated along the time.
- [std_serum_files.txt](data/std_serum_files.txt): table containing the mzML
  files from the re-measurement of the standard mixes in January 2020. This
  contains:
  - measurement of standard mixes in water.
  - measurement of standard mixes in serum (QC CHRIS Pool samples).
  All mixes are added in a *low* and *high* concentration, for the high
  concentrations also DDA MS2 data was measured with two different collision
  energies.


## Contributors

- Marilyn De Graeve
- Andrea Vicini
- Vinicius Verri Hernandes
- Mar Garcia-Aloy
- Johannes Rainer (contact)


## *Old* R files

- [preprocessing-standards.Rmd](preprocessing-standards.Rmd): preprocessing of
  LC-MS(/MS) data for the standards data set.
- [EIC_superposed.Rmd](EIC_superposed.Rmd): Markdown file for printing 
  the EIC of one compound using all samples. The plots are saved 
  in the folder `plots`.
- [dinamic_range.Rmd](dinamic_range.Rmd): Markdown file for getting the 
  calibration curves (and their corresponding plots) of each standard in 
  each ionization mode (i.e., polarity). Plots are saved in the folder 
  `images/linearity` and data in the text file `dinamic_range_POL.txt`.
- [get_rt.Rmd](get_rt.Rmd): Markdown file for obtaining the RTs of the 
  maximum peak corresponding to each compound (considering its exact mass). 
  This serve as a starting point for define the RT of each compound. 
  At the begining we have to specify the following parameters:  
  - `study`: specify which of the 2 studies we want to focus 
  (i.e., "internal_standards" or "standards_dilution").
  - `mixnum`: in case we specified `study <- "standards_dilution"`, 
  specify to which MIX we want to focus.
  - `polarity`: specify in which polarity ("POS" or "NEG") we want to focus.
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
