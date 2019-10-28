NEG = c(
  # Formic Acid:
  44.9820, (46.00547*2)-1.007276, (44.9820 + seq(14)*67.98745),
  (44.9820+(c(3,4)*67.98745)+57.95865),
  
  # ??
  623.3915, (623.3915+seq(4)*67.98745)
  )

POS = c(124.0868, 90.97639, 90.97639+(seq(13)*67.98745))

max.len = max(length(NEG), length(POS))

NEG = c(NEG, rep(NA, max.len - length(NEG)))
POS = c(POS, rep(NA, max.len - length(POS)))

exclusion_mz <- data.frame(NEG, POS)

write.table(exclusion_mz, "data/exclusion_mz.txt", sep="\t", row.names = FALSE)
