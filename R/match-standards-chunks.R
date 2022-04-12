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
setMSnbaseFastLoad(FALSE)
register(SerialParam())
source("R/match-standards-functions.R")


## ---- general-settings ----
MIX_NAME <- paste0("Mix", ifelse(MIX < 10, paste0(0, MIX), MIX)) 
IMAGE_PATH <- paste0("images/match-standards-serum-", tolower(MIX_NAME),"/")
if (dir.exists(IMAGE_PATH)) unlink(IMAGE_PATH, recursive = TRUE)
RDATA_PATH <- paste0("data/RData/match-standards-serum-", tolower(MIX_NAME), "/")
dir.create(IMAGE_PATH, showWarnings = FALSE, recursive = TRUE)
dir.create(RDATA_PATH, showWarnings = FALSE, recursive = TRUE)
#' Define the mzML files *base* path (/data/massspec/mzML/ on the cluster)
MZML_PATH <- "~/mix01/"
#' MZML_PATH <- "/data/massspec/mzML/"
ALL_NL_MATCH <- FALSE                   # run maching agains neutral loss db
library(knitr)
opts_chunk$set(cached = FALSE, message = FALSE, warning = FALSE,
               fig.width = 10, fig.height = 8)
#'               dev = c("png", "pdf"), fig.path = IMAGE_PATH)


## ---- standards-table ----
std_dilution <- read.table("data/standards_dilution.txt",
                           sep = "\t", header = TRUE)
std_dilution$exactmass <- calculateMass(std_dilution$formula)
colnames(std_dilution)[colnames(std_dilution) == "HMDB.code"] <- "HMDB"
std_dilution01 <- std_dilution[std_dilution$mix == MIX, ]
pandoc.table(std_dilution01[, c("name", "formula", "HMDB", "RT", "POS", "NEG")], 
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
library(CompoundDb)

cdb <- CompDb("data/CompDb.Hsapiens.HMDB.5.0.sqlite")
hmdb_pos <- hmdb_neg <- Spectra(cdb)
hmdb_pos$precursorMz <- mass2mz(hmdb_pos$exactmass, adduct = "[M+H]+")[, 1L]
hmdb_neg$precursorMz <- mass2mz(hmdb_pos$exactmass, adduct = "[M-H]-")[, 1L]
#' Neutral loss spectra
hmdb_pos_nl <- neutralLoss(hmdb_pos, param = PrecursorMzParam())
hmdb_neg_nl <- neutralLoss(hmdb_neg, param = PrecursorMzParam())
#' HMDB with only MS2 for current standards
hmdb_std <- Spectra(cdb, filter = ~ compound_id == std_dilution01$HMDB)

library(MsBackendMassbank)
library(RSQLite)
con <- dbConnect(SQLite(), "data/MassBank.sqlite")
mbank <- Spectra(con, source = MsBackendMassbankSql())
mbank$name <- mbank$compound_name
#' Neutral loss spectra
mbank_nl <- mbank[!is.na(mbank$precursorMz)]
mbank_nl <- neutralLoss(mbank_nl, param = PrecursorMzParam())


## ---- all-ms2 ----
fls <- fls[grep("_CE", fls)]

all_ms2 <- Spectra(fls, source = MsBackendMzR()) |>
filterMsLevel(msLevel. = 2L) |>
filterIntensity(intensity = low_int)

all_ms2 <- all_ms2[lengths(all_ms2) > 1]
all_ms2 <- addProcessing(all_ms2, scale_int)
all_ms2 <- setBackend(all_ms2, MsBackendDataFrame())
all_ms2 <- applyProcessing(all_ms2)


## ---- compare-spectra-param
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
    tmp, std_dilution01[!std_dilution01$HMDB %in% tmp$HMDB, c("name", "HMDB")])
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
             "[M+2Na-H]+", "[M+2K-H]+")

## ---- define-adducts-neg ----
adducts <- adducts("negative")[c("[M-H]-", "[M+Cl]-"),
                               c("mass_multi", "mass_add")] 
adducts <- rbind(adducts, "[M+HCOO]-" = c(1, calculateMass("HCOO"))) 

## ---- match-features ----
prm <- Mass2MzParam(adducts = adducts, ppm = 30)
fmat <- c(featureDefinitions(data), ttest)
mtchs <- matchMz(fmat, std_dilution01, param = prm, mzColname = "mzmed")
mtchs_sub <- mtchs[whichQuery(mtchs)]
mD <- matchedData(
    mtchs_sub, columns = c("mzmed", "ppm_error", "rtmed", "target_RT",
                           "target_name", "target_HMDB", "adduct", 
                           "high_low_diff", "pvalue", "mean_high",
                           "mean_low"))
mD <- mD[-which(mD$high_low_diff < 1), ]
mD <- mD[!(is.na(mD$mean_high) & !is.na(mD$mean_low)), ]
mD <- mD[order(mD$target_name), ]


## ---- table-standard-no-feature ----
pandoc.table(std_dilution01[!std_dilution01$name %in% mD$target_name,
                            c("name", "HMDB", "formula", "RT", "POS", "NEG")],
             style = "rmarkdown", split.tables = Inf,
             caption = "Not matched standards")


## ---- table-feature-matches ----
tmp <- as.data.frame(mD[mD$target_name == std, ])
tmp <- tmp[order(tmp$rtmed), ]
tmp$feature_group <- group_features(data_FS, rownames(tmp))
std_ms2 <- extract_ms2(data, rownames(tmp))
tmp$n_ms2 <- 0
if(length(std_ms2)) {
   nms2 <- table(std_ms2$feature_id)
   tmp[names(nms2), "n_ms2"] <- unname(nms2) 
}
pandoc.table(tmp[, c("mzmed", "ppm_error", "rtmed", "target_RT",
                     "feature_group", "adduct", "n_ms2", 
                     "mean_high", "mean_low")], style = "rmarkdown",
             split.tables = Inf,
             caption = paste0("Feature to standards matches for ", std))


## ---- plot-ms2-hmdb-heatmap ----
hmdb_id <- std_dilution$HMDB[std_dilution$name == std]
std_hmdb <- Spectra(cdb, filter = ~ compound_id == hmdb_id)
if (length(std_hmdb)) {
    sim <- compareSpectra(std_ms2, std_hmdb, ppm = csp@ppm,
                          tolerance = csp@tolerance)
    ann <- data.frame(feature_id = std_ms2$feature_id, rt = rtime(std_ms2))
    if (is.matrix(sim) && all(dim(sim) > 1)) {
        rownames(ann) <- rownames(sim)
        pheatmap(sim, annotation_row = ann, breaks = seq(0, 1, length.out = 101),
                 color = colorRampPalette((brewer.pal(n = 7, name = "YlOrRd")))(100))
    } else cat("Similarities: ", sim)
} else cat("No reference spectrum in HMDB")

## ---- plot-ms2-mbank-heatmap ----
std_mbank <- get_mbank(mbank, inchikey = unique(std_hmdb$inchikey))
if (length(std_mbank)) {
    sim <- compareSpectra(std_ms2, std_mbank, ppm = csp@ppm,
                          tolerance = csp@tolerance)
    ann <- data.frame(feature_id = std_ms2$feature_id, rt = rtime(std_ms2))
    if (is.matrix(sim) && all(dim(sim) > 1)) {
        rownames(ann) <- rownames(sim)
        pheatmap(sim, annotation_row = ann, breaks = seq(0, 1, length.out = 101),
                 color = colorRampPalette((brewer.pal(n = 7, name = "YlOrRd")))(100))
    } else cat("Similarities: ", sim)
} else cat("No reference spectrum in MassBank")

## ---- plot-ms2-hmdb-nl-heatmap ----
std_hmdb_nl <- hmdb_nl[hmdb_nl$compound_id == hmdb_id]
std_ms2_nl <- neutralLoss(std_ms2, PrecursorMzParam())
sim <- compareSpectra(std_ms2_nl, std_hmdb_nl, ppm = csp_nl@ppm,
                      tolerance = csp_nl@tolerance)
ann <- data.frame(feature_id = std_ms2_nl$feature_id, rt = rtime(std_ms2_nl))
if (is.matrix(sim) && all(dim(sim) > 1)) {
rownames(ann) <- rownames(sim)
pheatmap(sim, annotation_row = ann, breaks = seq(0, 1, length.out = 101),
         color = colorRampPalette((brewer.pal(n = 7, name = "YlOrRd")))(100))
} else cat("Similarities: ", sim)

## ---- plot-ms2-mbank-nl-heatmap ----
std_mbank_nl <- get_mbank(mbank, inchikey = unique(std_hmdb$inchikey), nl = TRUE)
std_ms2_nl <- neutralLoss(std_ms2, PrecursorMzParam())
sim <- compareSpectra(std_ms2_nl, std_mbank_nl, ppm = csp_nl@ppm,
                      tolerance = csp_nl@tolerance)
ann <- data.frame(feature_id = std_ms2$feature_id, rt = rtime(std_ms2))
if (is.matrix(sim) && all(dim(sim) > 1)) {
rownames(ann) <- rownames(sim)
pheatmap(sim, annotation_row = ann, breaks = seq(0, 1, length.out = 101),
         color = colorRampPalette((brewer.pal(n = 7, name = "YlOrRd")))(100))
} else cat("Similarities: ", sim)
