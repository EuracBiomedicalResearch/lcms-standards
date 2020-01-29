neg_FA <- c(
  46.00547-1.007276,     # [M-H]-
  (46.00547*2)-1.007276, # [2M-H]-
  (46.00547-1.007276 + seq(14)*67.98745), # adducts of CHOONa (sodium formate)
  (46.00547-1.007276 + (c(3,4)*67.98745) + 57.95865), # adducts of CHOONa & NaCl (sodium chloride)
  758.86518, 758.86518 + (seq(3)*67.98745)) # mz 758.86 (?) + adducts of CHOONa

neg_IS_030 <- c(61.98875, 115.9206, 116.92741, 132.9233, 144.9235, 
                383.18958, 384.18454, 457.22673, 540.2722, 540.2722+1.00335)
neg_01_030 <- c(248.97233)
neg_IS_110 <- c(277.18426)
neg_IS_122 <- c(311.16875, 325.18432, 339.19983)
neg_IS_178 <- c(195.8112, 197.8083, 199.80514)
neg_XX <- c(623.3915, (623.3915+seq(4)*67.98745))

NEG_ions <- c(neg_FA, neg_IS_030, neg_01_030, 
              neg_IS_110, neg_IS_122, neg_IS_178, neg_XX)
NEG_cmts <- c(rep("Formic Acid",    length(neg_FA)),
              rep("IS - RT 030",    length(neg_IS_030)),
              rep("MIX01 - RT 030", length(neg_01_030)),
              rep("IS - RT 110",    length(neg_IS_110)),
              rep("IS - RT 122",    length(neg_IS_122)),
              rep("IS - RT 178",    length(neg_IS_178)),
              rep("XX",             length(neg_XX)))

##############################################################################

pos_XX <- c(
  124.0868, # mz 124.09 (?)
  90.97639, 90.97639 + (seq(13)*67.98745)) # mz 90.97 (?) + adducts of CHOONa
pos_IS_026 <- c(223.0595, 283.0262)
pos_IS_027 <- c(581.2343, 637.2957, (637.2957+1.00335), (637.2957+21.98198))
pos_IS_030 <- c(341.2999)
pos_IS_031 <- c(193.1589, 437.1869)
pos_IS_033 <- c(455.3272)
pos_01_040 <- c(60.04581)
pos_01_045 <- c(74.06112)
pos_01_074 <- c(114.0662, (114.0662+21.98198))
pos_IS_121 <- c(135.09299)
pos_IS_196 <- c(90.50735, 184.9847)
pos_IS_199 <- c(89.50731, 113.96138, 141.95797)
pos_IS_203 <- c(182.9841)

POS_ions <- c(pos_XX, pos_IS_026, pos_IS_027, pos_IS_030, pos_IS_031, 
              pos_IS_033, pos_01_040, pos_01_045, pos_01_074, pos_IS_121, 
              pos_IS_196, pos_IS_199, pos_IS_203)

POS_cmts <- c(rep("XX",          length(pos_XX)),
              rep("IS - RT 026", length(pos_IS_026)),
              rep("IS - RT 027", length(pos_IS_027)),
              rep("IS - RT 030", length(pos_IS_030)),
              rep("IS - RT 031", length(pos_IS_031)),
              rep("IS - RT 033", length(pos_IS_033)),
              rep("IS - RT 040", length(pos_01_040)),
              rep("IS - RT 045", length(pos_01_045)),
              rep("IS - RT 074", length(pos_01_074)),
              rep("IS - RT 121", length(pos_IS_121)),
              rep("IS - RT 196", length(pos_IS_196)),
              rep("IS - RT 199", length(pos_IS_199)),
              rep("IS - RT 203", length(pos_IS_203)))

##############################################################################

max.len = max(length(NEG_ions), length(POS_ions))

NEG_ions = c(NEG_ions, rep(NA, max.len - length(NEG_ions)))
NEG_cmts = c(NEG_cmts, rep(NA, max.len - length(NEG_cmts)))
POS_ions = c(POS_ions, rep(NA, max.len - length(POS_ions)))
POS_cmts = c(POS_cmts, rep(NA, max.len - length(POS_cmts)))

exclusion_mz <- data.frame(NEG_ions, NEG_cmts, POS_ions, POS_cmts)

write.table(exclusion_mz, "data/exclusion_mz.txt", sep="\t", row.names = FALSE)
