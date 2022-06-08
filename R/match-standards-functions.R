#' Functions to process the spectra
low_int <- function(x, ...) {
    x > max(x, na.rm = TRUE) * 0.05
}
scale_int <- function(x, ...) {
    maxint <- max(x[, "intensity"], na.rm = TRUE)
    x[, "intensity"] <- 100 * x[, "intensity"] / maxint
    x
}

#' helper function to group features and return their feature group IDs.
group_features <- function(x, features, groupEic = FALSE) {
    featureGroups(x) <- rep(NA, nrow(featureDefinitions(x)))
    featureGroups(x)[rownames(featureDefinitions(x)) %in% features] <- "FG"
    x <- groupFeatures(x, param = SimilarRtimeParam(5))
    x <- groupFeatures(x, param = AbundanceSimilarityParam(threshold = 0.6,
                                                           transform = log2))
    if (groupEic)
        x <- groupFeatures(x, param = EicSimilarityParam(threshold = 0.5,
                                                         n = 3))
    featureGroups(x)[match(features, rownames(featureDefinitions(x)))]
}

#' Plot EICs.
plot_eics <- function(x, std, tab, MP, std_ms2) {
    eics <- featureChromatograms(x, features = rownames(tab))
    dr <- file.path(IMAGE_PATH, std)
    dir.create(dr, showWarnings = FALSE)
    col_sample <- col_group[eics$group]

    for (i in seq_len(nrow(eics))) {
        ft <- rownames(tab)[i]
        fl <- file.path(dr, paste0(MP,"_", ft, ".png"))
        png(fl, width = 12, height = 8, units = "cm", res = 300, pointsize = 6)
        eic <- eics[i, ]
        plot(eic, col = paste0(col_sample, 80),
             peakBg = paste0(col_sample[chromPeaks(eic)[, "sample"]], 40))
        abline(v = tab$target_RT[1], lty = 2, col = "#00000060")
        grid()
        legend("topleft", c(std, tab$adduct[i], ft))
        legend("topright", col = col_group, legend = names(col_group), lty = 1)
        if (length(std_ms2) && length(idx <- which(std_ms2$feature_id == ft)))
            abline(v = rtime(std_ms2)[idx], col = col_group["MSMS"], lty = 2)
        dev.off()
    }
}

plot_spectra <- function(x) {
    if (length(x)) {
        xl <- range(unlist(mz(x)))
        lx <- length(x)
        par(mfrow = c(ceiling(sqrt(lx)), round(sqrt(lx))),
            mar = c(4.2, 4.5, 1.5, 0.5))
        for (i in seq_len(lx)) {
            plotSpectra(x[i], xlim = xl)
            legend("top", bg = NA, horiz = TRUE,
                   legend = c(x$feature_id[i],
                              spectraNames(x)[i]))
        }
    } else cat("No MS2 spectra.")
}

extract_ms2 <- function(x, features, ppm = 20, tolerance = 0) {
    res <- featureSpectra(x, expandRt = 3, return.type = "Spectra",
                          features = features, ppm = ppm, expandMz = tolerance)
    res <- filterIntensity(res, intensity = low_int)
    res <- res[lengths(res) > 1]
    if (!length(res))
        return(res)
    res <- addProcessing(res, scale_int)
    spectraNames(res) <- paste0(res$feature_id, "_", spectraNames(res))
    res
}

#' helper function that selects experimental MS2 spectra (std_ms2) with a
#' similarity higher than the specified cutoff, creates mirror plots and
#' returns the selected `Spectra`.
#'
#' This function uses *global* variables `sim`, `std_ms2` and `hmdb`
plot_select_ms2 <- function(query, target, cutoff,
                            ppm = 40, tolerance = 0) {
    if (is.null(dim(sim)))
        idx <- which(matrix(sim, length(query), length(target)) > cutoff,
                     arr.ind = TRUE)
    else idx <- which(sim > cutoff, arr.ind = TRUE)
    if (nrow(idx)) {
        par(mfrow = c(round(sqrt(nrow(idx))), ceiling(sqrt(nrow(idx)))))
        for (i in seq_len(nrow(idx))) {
            a <- idx[i, 1L]
            b <- idx[i, 2L]
            plotSpectraMirror(query[a], addProcessing(target[b], scale_int),
                              main = paste(spectraNames(query)[a], target$name[b]),
                              ppm = ppm, tolerance = tolerance)
        }
        query[unique(idx[, "row"])]
    } else cat("No spectra with selected similarity")
}

match_table <- function(x, sv = c("rtime", "target_name",
                                  "score")) {
    x <- x[whichQuery(x)]
    tab <- spectraData(x, c("feature_id", sv))
    tab$precursorMz_diff <- x$precursorMz - x$target_precursorMz
    rn <- strsplit(rownames(tab), "_")
    rn <- vapply(rn, function(z) z[2L], character(1))
    tab$spectrum_id <- rn
    fids <- split(tab$feature_id, tab$spectrum_id)
    fids <- vapply(fids, function(z) paste0(unique(z), collapse = ";"),
                   character(1))
    tab <- as.data.frame(tab)
    rownames(tab) <- NULL
    tab <- unique(tab[, c(sv, "precursorMz_diff", "spectrum_id")])
    tab$feature_id <- fids[tab$spectrum_id]
    tab
}

summarize_match_table <- function(x, name = "target_name") {
  f <- paste(x[, name], x$spectrum_id)
  x <- lapply(split(x, f), function(z)
    z[order(z$score, decreasing = TRUE), ][1, ])
  x <- do.call(rbind, x)
  x <- x[order(x$score, decreasing = TRUE), ]
  x$selected <- rep("", nrow(x))
  if(length(std_ms2_sel)) {
    sel_ids <- vapply(strsplit(spectraNames(std_ms2_sel), "_"),
                      function(z) z[2], character(1))
    x$selected[x$spectrum_id %in% sel_ids] <- "X"
  }
  rownames(x) <- NULL
  x
}

perform_match <- function(query, target, sv, name, param, similarity = 0.7) {
    mtch <- matchSpectra(
        query, target, param = param)
    if (length(whichQuery(mtch))) {
        st <- match_table(mtch, sv) |>
        summarize_match_table(name)
        pandoc.table(
            st[st$score > similarity, ], style = "rmarkdown", split.table = Inf,
            caption = paste0("Experimental MS2 spectra with best match to MS2",
                             " spectra of different compounds. The best match ",
                             "per compound is reported. Rows are ordered by ",
                             "similarity. Only matches with a similarity ",
                             "larger than ", similarity, " are reported."))
    } else cat("No matches found.")
}

#' get mbank spectra for a standard identified with an inchikey.
get_mbank <- function(x, inchikey = character(), nl = FALSE) {
    x <- x[which(x$inchikey %in% inchikey)]
    if (length(x) && nl) {
        no_prec <- is.na(precursorMz(x))
        if (any(no_prec)) {
            idx <- which(no_prec & !is.na(x$adduct))
            x$precursorMz[idx] <- mass2mz(x$exactmass[idx],
                                          adduct = x$adduct[idx])
            x <- x[!is.na(precursorMz(x))]
        }
        x <- neutralLoss(x, param = PrecursorMzParam())
    }
    x
}

#' Calculates similarity between experimental MS2 spectra and reference
#' MS2 and returns the similarity matrix.
plot_ms2_similarity_heatmap <- function(x, y, param) {
    if (length(y)) {
        sim <- compareSpectra(x, y, ppm = param@ppm,
                              tolerance = param@tolerance,
                              simplify = FALSE)
        ann <- data.frame(feature_id = x$feature_id, rt = rtime(x))
        if (is.matrix(sim) && all(dim(sim) > 1)) {
            rownames(ann) <- rownames(sim)
            pheatmap(sim, annotation_row = ann,
                     breaks = seq(0, 1, length.out = 101),
                     color = colorRampPalette(
                     (brewer.pal(n = 7, name = "YlOrRd")))(100))
        } else cat("Similarities: ", sim)
        sim
    } else {
        cat("No reference spectra available")
        matrix()
    }
}
