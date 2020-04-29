library(ggplot2)
library(tidyr)
library(RColorBrewer)
load("tmp.RData")
tmp <- data.frame(formula_finder)
tmp$abbreviation <- rownames(tmp)
tmp <- merge(std_info, tmp, by = "abbreviation")
tmp <- tmp[tmp$abbreviation != "NADH", ]
tmp$mzgroup <- 1
tmp$mzgroup[tmp$mzneut > 250] <- 2
tmp$mzgroup[tmp$mzneut > 600] <- 3

tmplong <- gather(tmp,
                   key = "property",
                   value = "value",
                   start, thr_filters, A1_filter, C1_filter, A2_filter, rank, rankp)
tmplong$value[is.na(tmplong$value)] <- 0
tmplong$property <- 
  factor(tmplong$property, 
         levels = c("start", "thr_filters", 
                    "A1_filter", "C1_filter", "A2_filter", "rank", "rankp"))
ggplot(tmplong, aes(x=property,y=value, group = abbreviation)) + 
  geom_point()  + geom_line(aes(col = brewer.pal(3, "Set1")[mzgroup])) + 
  theme(legend.position = "none")

tmplong <- gather(tmp,
                  key = "property",
                  value = "value",
                   A2_filter, rank, rankp)
tmplong$value[is.na(tmplong$value)] <- 0
tmplong$property <- 
  factor(tmplong$property, 
         levels = c("A2_filter", "rank", "rankp"))
ggplot(tmplong, aes(x=property,y=value, group = abbreviation)) + 
  geom_point()  + geom_line(aes(col = brewer.pal(3, "Set1")[mzgroup])) + 
  theme(legend.position = "none")
ggplot(tmplong, aes(x=property,y=value, group = abbreviation)) + 
  geom_point()  + geom_line(aes(col = brewer.pal(3, "Set1")[mzgroup])) + 
  theme(legend.position = "none") + ylim(0,200)
ggplot(tmplong, aes(x=property,y=value, group = abbreviation)) + 
  geom_point()  + geom_line(aes(col = brewer.pal(3, "Set1")[mzgroup])) + 
  theme(legend.position = "none") + ylim(0,70)
ggplot(tmplong, aes(x=property,y=value, group = abbreviation)) + 
  geom_point()  + geom_line(aes(col = brewer.pal(3, "Set1")[mzgroup])) + 
  theme(legend.position = "none") + ylim(0,10)

tmp$startp <- 100
tmp$thr_filtersp <- (tmp$thr_filters*100)/tmp$start
tmp$A1_filterp <- (tmp$A1_filter*100)/tmp$start
tmp$C1_filterp <- (tmp$C1_filter*100)/tmp$start
tmp$A2_filterp <- (tmp$A2_filter*100)/tmp$start
tmplong <- gather(tmp,
                  key = "property",
                  value = "value",
                  startp, thr_filtersp, A1_filterp, C1_filterp, A2_filterp)
tmplong$value[is.na(tmplong$value)] <- 0
tmplong$property <- 
  factor(tmplong$property, 
         levels = c("startp", "thr_filtersp", 
                    "A1_filterp", "C1_filterp", "A2_filterp"))
ggplot(tmplong, aes(x=property,y=value, group = abbreviation)) + 
  geom_point()  + geom_line(aes(col = brewer.pal(3, "Set1")[mzgroup])) + 
  theme(legend.position = "none")
tmplong <- tmplong[tmplong$value < 100, ]




tmpx <- tmpx[order(tmpx$A),]
tmpx$rank <- seq(nrow(tmpx))
tmpx <- tmpx[order(tmpx$Ap),]
tmpx$rankp <- seq(nrow(tmpx))
