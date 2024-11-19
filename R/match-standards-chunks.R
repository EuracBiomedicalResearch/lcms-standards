#' These are code blocks (chunks) used across the various match-standards Rmd
#' files. Defining them once here and re-using them across files helps avoid copy-paste errors

## ---- libraries ----
library(xcms)
library(pander)
library(MetaboCoreUtils)
library(MsCoreUtils)
library(MetaboAnnotation)
library(BiocParallel)
library(Spectra)
library(RColorBrewer)
library(MetaboAnnotation)
library(pheatmap)
library(MsFeatures)
library(MSnbase)
# setMSnbaseFastLoad(TRUE)
setMSnbaseFastLoad(FALSE)
register(SerialParam())
source("R/match-standards-functions.R")


## ---- general-settings ----
MIX_NAME <- paste0("Mix", ifelse(MIX < 10, paste0(0, MIX), MIX)) 
IMAGE_PATH <- paste0("images/match-standards-", tolower(MATRIX), "-",
                     tolower(MIX_NAME), "/")
#if (dir.exists(IMAGE_PATH)) unlink(IMAGE_PATH, recursive = TRUE)
RDATA_PATH <- paste0("data/RData/match-standards-", tolower(MATRIX), "-",
                     tolower(MIX_NAME), "/")
dir.create(IMAGE_PATH, showWarnings = FALSE, recursive = TRUE)
dir.create(RDATA_PATH, showWarnings = FALSE, recursive = TRUE)
#' Define the mzML files *base* path (/data/massspec/mzML/ on the cluster)
MZML_PATH <- "/Volumes/PortableSSD/mzML/"
MZML_PATH <- "/data/massspec/mzML/"
ALL_NL_MATCH <- FALSE                   # run matching against neutral loss db
library(knitr)
opts_chunk$set(cached = FALSE, message = FALSE, warning = FALSE,
               fig.width = 10, fig.height = 8)
#'               dev = c("png", "pdf"), fig.path = IMAGE_PATH)


## ---- standards-table ----
std_dilution <- read.table("data/standards_dilution.txt",
                           sep = "\t", header = TRUE)
std_dilution$exactmass <- calculateMass(std_dilution$formula)
colnames(std_dilution)[colnames(std_dilution) == "HMDB.code"] <- "HMDB"
std_dilution <- std_dilution[std_dilution$mix == MIX, ]
pandoc.table(std_dilution[, c("name", "formula", "HMDB", "RT", "POS", "NEG")], 
             style = "rmarkdown", split.tables = Inf,
             caption = paste0("Standards of ", MIX_NAME, ". Columns RT, POS ",
                              "and NEG contain expected retention times and ",
                              "adducts from a previous manual analysis."))

## ---- data-import ----
std_files <- read.table("data/std_serum_files.txt", header = TRUE)
std_files$concentration <- "blank"
std_files$concentration[grep("High", std_files$class)] <- "high"
std_files$concentration[grep("Low", std_files$class)] <- "low"
std_files$group <- std_files$concentration
std_files$group[grep("CE", std_files$mode)] <- "MSMS"

#' Subset and load data.
std_files <- std_files[which(std_files$type == MIX_NAME), ]
std_files$matrix <- "Water"
std_files$matrix[grep("^QC", std_files$class)] <- "Serum"
std_files <- std_files[std_files$matrix == MATRIX, ]
fls <- paste0(MZML_PATH, std_files$folder, "/", std_files$mzML)
data_all <- readMSData(fls, pdata = new("NAnnotatedDataFrame", std_files),
                       mode = "onDisk")
data_all <- filterRt(data_all, rt = c(0, 350))
data_all <- filterEmptySpectra(data_all)
#' Define colors
col_group <- brewer.pal(9, "Set1")[c(1, 5, 4, 9)]
names(col_group) <- c("high", "low", "MSMS", "blank")
col_group <- col_group[unique(data_all$group)]


## ---- load-reference-databases ----
#' Maybe we should/need to clean also the reference database the same way we
#' do with the experimental MS2 spectra, i.e. clean them to remove low
#' intensity peaks.
library(CompoundDb)

cdb <- CompDb("data/CompDb.Hsapiens.HMDB.5.0.sqlite")
hmdb_pos <- hmdb_neg <- Spectra(cdb)
hmdb_pos$precursorMz <- mass2mz(hmdb_pos$exactmass, adduct = "[M+H]+")[, 1L]
hmdb_neg$precursorMz <- mass2mz(hmdb_pos$exactmass, adduct = "[M-H]-")[, 1L]
#' Neutral loss spectra
nl_param <- PrecursorMzParam(filterPeaks = "removePrecursor", ppm = 50,
                             tolerance = 0)
hmdb_pos_nl <- neutralLoss(hmdb_pos, param = nl_param)
hmdb_pos_nl <- hmdb_pos_nl[lengths(hmdb_pos_nl) > 0]
hmdb_neg_nl <- neutralLoss(hmdb_neg, param = nl_param)
hmdb_neg_nl <- hmdb_neg_nl[lengths(hmdb_neg_nl) > 0]
#' HMDB with only MS2 for current standards
hmdb_std <- Spectra(cdb, filter = ~ compound_id == std_dilution$HMDB)

## library(MsBackendMassbank)
## library(RSQLite)
#con <- dbConnect(SQLite(), "data/MassBank.sqlite") #NO, alternatively do:
#con <- dbConnect(SQLite(), "data/MassBank.sql") #error
library(AnnotationHub)
ah <- AnnotationHub()
con <- ah[["AH116166"]]

mbank <- Spectra(con)
#mbank$name <- mbank$compound_name   #already present in AH116166
#' Neutral loss spectra
mbank_nl <- mbank[!is.na(mbank$precursorMz)]
mbank_nl <- neutralLoss(mbank_nl, param = nl_param)
mbank_nl <- mbank_nl[lengths(mbank_nl) > 0]

## ---- setup-ion-db ----
fl <- "data/IfB_HILIC.IonDb.HMDB.5.0.sqlite"
if (!file.exists(fl)) {
    idb <- IonDb(fl, cdb)
} else {
    ## Delete any ions for the current mix/matrix
    idb <- IonDb(fl)
    all_ions <- ions(idb, columns = ionVariables(idb, includeId = TRUE))
    ids <- all_ions$ion_id[all_ions$sample_matrix == MATRIX &
                           all_ions$original_sample == paste0("Std_", MIX_NAME)]
    if (length(ids)) idb <- deleteIon(idb, ids = ids)
    ## Delete any MS2 spectra from the current files
    all_sps <- Spectra(idb)
    ids <- all_sps$spectrum_id[all_sps$original_file %in% std_files$mzML]
    if (length(ids)) idb <- deleteSpectra(idb, ids = ids)
}

## ---- all-ms2 ----
fls <- fls[grep("_CE", fls)]

all_ms2 <- Spectra(fls, source = MsBackendMzR()) |>
filterMsLevel(msLevel. = 2L) |>
filterIntensity(intensity = low_int)

all_ms2 <- all_ms2[lengths(all_ms2) > 1]
all_ms2 <- addProcessing(all_ms2, scale_int)
all_ms2 <- setBackend(all_ms2, MsBackendDataFrame())
all_ms2 <- applyProcessing(all_ms2)


## ---- compare-spectra-param ----
## Settings for matchSpectra
csp <- CompareSpectraParam(ppm = 40, requirePrecursor = FALSE)
## Settings for matchSpectra with neutral loss spectra
csp_nl <- CompareSpectraParam(ppm = 0, tolerance = 0.1,
                              requirePrecursor = FALSE)


## ---- table-all-ms2-hmdb ----
all_match <- matchSpectra(
    hmdb_std, all_ms2,
    param = csp)
all_match <- all_match[whichQuery(all_match)]

tab <- matchedData(all_match, c("name", "compound_id", "score", "target_rtime"))

tmp <- split(as.data.frame(tab), tab$compound_id)
tmp <- do.call(rbind, lapply(tmp, function(z) {
    data.frame(name = z$name[1L],
               HMDB = z$compound_id[1L],
               ms2 = nrow(z),
               mean_rt = mean(z$target_rtime),
               min_rt = min(z$target_rtime),
               max_rt = max(z$target_rtime),
               mean_sim = mean(z$score))
}))
tmp <- rbindFill(
    tmp, std_dilution[!std_dilution$HMDB %in% tmp$HMDB, c("name", "HMDB")])
rownames(tmp) <- NULL
pandoc.table(tmp[order(tmp$name), ], style = "rmarkdown", split.table = Inf,
             caption = "Standards with matching reference MS2 spectra.")


## ---- preprocessing ----
cwp <- CentWaveParam(ppm = 50,
                     peakwidth = c(2, 18),
                     snthresh = 5,
                     mzdiff = 0.001,
                     prefilter = c(4, 300),
                     noise = 100,
                     integrate = 2)
data <- findChromPeaks(data, param = cwp) 
#' Peak refinement
mnp <- MergeNeighboringPeaksParam(expandRt = 3.5, expandMz = 0.001,
                                  minProp = 3/4)
data <- refineChromPeaks(data, param = mnp)
#' Alignment
pdp1 <- PeakDensityParam(sampleGroups = data$matrix, bw = 3,
                         minFraction = 0.7, binSize = 0.015)
data <- groupChromPeaks(data, param = pdp1)
pgp <- PeakGroupsParam(minFraction = 0.8, extraPeaks = 100, span = 0.8)
data <- adjustRtime(data, param = pgp)
#' Correspondence analysis
pdp2 <- PeakDensityParam(sampleGroups = data$mode, bw = 3,
                         minFraction = 0.3, binSize = 0.015)
data <- groupChromPeaks(data, param = pdp2)
#' Gap-filling
data <- fillChromPeaks(data, param = ChromPeakAreaParam())
save(data, file = paste0(RDATA_PATH, "processed_data_", POLARITY, ".RData"))


## ---- abundance-difference ----
fVlog2 <- log2(featureValues(data, value = "into",
                             method = "sum", filled = TRUE))
high <- grep("High", colnames(fVlog2))
high <- high[-grep("CE", colnames(fVlog2)[high])]
low <- grep("Low", colnames(fVlog2))
ttest <- t(apply(fVlog2, 1, function(x) {
    a <- x[high]
    b <- x[low]
    if (sum(!is.na(a)) > 1 && sum(!is.na(b)) > 1) {
        res <- t.test(a, b, mu = 0)
        c(high_low_diff = unname(res$estimate[1] - res$estimate[2]),
          pvalue = res$p.value, mean_high = mean(a, na.rm = TRUE),
          mean_low = mean(b, na.rm = TRUE))
    } else c(high_low_diff = mean(a, na.rm = TRUE) - mean(b, na.rm = TRUE),
             pvalue = NA_real_, mean_high = mean(a, na.rm = TRUE),
             mean_low = mean(b, na.rm = TRUE))
}))


## ---- define-adducts-pos ----
adducts <- c("[M+H]+", "[M+2H]2+", "[M+Na]+", "[M+K]+", "[M+NH4]+",
             "[M+H-H2O]+", "[M+H+Na]2+", "[M+2Na]2+", "[M+H-NH3]+",
             "[M+2Na-H]+", "[M+2K-H]+",
             "[2M+H]+", "[M+H-H4O2]+", "[M+H-Hexose-H2O]+", "[M+H-CH2O2]+")

## ---- define-adducts-neg ----
adducts <- adducts("negative")[c("[M-H]-", "[M+Cl]-", "[M-H+HCOONa]-",
                                 "[2M-H]-", "[M+CHO2]-"),
                               c("mass_multi", "mass_add")] 

## ---- match-features ----
prm <- Mass2MzParam(adducts = adducts, ppm = 30)
fmat <- c(featureDefinitions(data), ttest)
mtchs <- matchMz(fmat, std_dilution, param = prm, mzColname = "mzmed")
mtchs_sub <- mtchs[whichQuery(mtchs)]
mD <- matchedData(
    mtchs_sub, columns = c("mzmed", "ppm_error", "rtmed", "target_RT",
                           "target_name", "target_HMDB", "adduct", 
                           "high_low_diff", "pvalue", "mean_high",
                           "mean_low"))
mD <- mD[-which(mD$high_low_diff < 0.7), ]
mD <- mD[which(abs(mD$rtmed - mD$target_RT) < 10), ]  #max 10s rt deviation, 2-sided
mD <- mD[order(mD$target_name), ]


## ---- table-standard-no-feature ----
pandoc.table(std_dilution[!std_dilution$name %in% mD$target_name,
                            c("name", "HMDB", "formula", "RT", "POS", "NEG")],
             style = "rmarkdown", split.tables = Inf,
             caption = "Not matched standards")


## ---- table-feature-matches ----
feature_table <- as.data.frame(mD[mD$target_name == std, ])
feature_table <- feature_table[order(feature_table$rtmed), ]
feature_table$feature_group <- group_features(data_FS, rownames(feature_table))
std_ms2 <- extract_ms2(data, rownames(feature_table), ppm = FEATURE_MS2_PPM,
                       tolerance = FEATURE_MS2_TOLERANCE)
feature_table$n_ms2 <- 0
if(length(std_ms2)) {
   nms2 <- table(std_ms2$feature_id)
   feature_table[names(nms2), "n_ms2"] <- unname(nms2) 
}
pandoc.table(feature_table[, c("mzmed", "ppm_error", "rtmed", "target_RT",
                               "feature_group", "adduct", "n_ms2", 
                               "mean_high", "mean_low")], style = "rmarkdown",
             split.tables = Inf,
             caption = paste0("Feature to standards matches for ", std))


## ---- add-ions ----
stopifnot(all(fts$feature_id %in% rownames(feature_table)))
fts$sample_matrix <- MATRIX
fts$original_sample <- paste0("Std_", MIX_NAME)
fts$ion_mz <- feature_table[fts$feature_id, "mzmed"]
fts$ion_rt <- round(feature_table[fts$feature_id, "rtmed"])
fts$ion_relative_intensity <- feature_table[fts$feature_id, "mean_high"]
fts$ion_relative_intensity <- fts$ion_relative_intensity /
    max(fts$ion_relative_intensity)
fts$ion_adduct <- feature_table[fts$feature_id, "adduct"]
fts$compound_id <- hmdb_id
fts$polarity <- POLARITY
idb <- insertIon(idb, fts, addColumns = TRUE)
pandoc.table(fts, split.tables = Inf, style = "rmarkdown")
rm(fts)

## ---- add-ms2-spectra ----
ms2$original_spectrum_id <- spectraNames(ms2)
ms2$compound_id <- hmdb_id
ms2$original_file <- basename(ms2$dataOrigin)
ms2$collisionEnergy <- 20
ms2$collisionEnergy[grep("CE30", ms2$original_file)] <- 30
ms2$predicted <- FALSE
ms2$instrument_type <- "LC-ESI-QTOF"
ms2$instrument <- "Sciex TripleTOF 5600+"
ms2$adduct <- feature_table[ms2$feature_id, "adduct"]
idb <- insertSpectra(
    idb, ms2, columns = c("original_spectrum_id", "compound_id", "polarity",
                          "collisionEnergy", "predicted", "instrument_type",
                          "instrument", "precursorMz", "adduct",
                          "original_file", "confidence", "acquisitionNum",
                          "rtime"))
spectraData(ms2, c("rtime", "original_file", "adduct", "confidence")) |>
    as.data.frame() |>
    pandoc.table(style = "rmarkdown", split.tables = Inf)
rm(ms2)

## ---- iondb-summary ----
fts <- ions(idb)
fts <- fts[fts$original_sample == paste0("Std_", MIX_NAME) &
           fts$sample_matrix == MATRIX, ]

fts <- split(fts, fts$compound_id)
fts <- do.call(rbind, lapply(fts, function(z) {
    data.frame(adducts = paste(z$ion_adduct, collapse = ", "),
               rt = mean(z$ion_rt))
}))
sps <- Spectra(idb)
sps <- sps[sps$original_file %in% basename(fileNames(data_all))]
n_ms2 <- table(sps$compound_id)
idx <- match(rownames(fts), std_dilution$HMDB)
fts <- cbind(name = std_dilution[idx, "name"],
             fts, n_ms2 = as.integer(n_ms2[rownames(fts)]),
             old_rt = std_dilution[idx, "RT"])
pandoc.table(fts[, c("name", "rt", "old_rt", "n_ms2", "adducts")],
             style = "rmarkdown", split.table = Inf)
