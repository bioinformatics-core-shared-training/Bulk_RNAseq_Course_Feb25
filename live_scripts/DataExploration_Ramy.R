library(tximport)
library(DESeq2)
library(tidyverse)

sampleinfo <- read_tsv("data/samplesheet.tsv", col_types = c("cccc"))
arrange(sampleinfo, Status, TimePoint, Replicate)


files <- file.path("salmon", sampleinfo$SampleName, "quant.sf")
files

files <- set_names(files, sampleinfo$SampleName)
files 

tx2gene <- read_tsv("references/tx2gene.tsv")

txi <- tximport(files, type= 'salmon', tx2gene = tx2gene)
str(txi)

head(txi$counts)

dir.create("salmon_outputs")
saveRDS(txi, file = "salmon_outputs/txi.rds")

rawCounts <- round(txi$counts, 0)
head(rawCounts)


dim(rawCounts)
keep <- rowSums(rawCounts) > 5

table(keep, useNA = "always")

filtCounts <- rawCounts[keep,]
dim(filtCounts)

summary(filtCounts)
boxplot(filtCounts, main = 'Raw counts', las = 2)

plot(rowMeans(filtCounts), rowSds(filtCounts),
     main = 'Raw counts: sd vs mean',
     xlim = c(0, 10000),
     ylim = c(0, 5000))

logCounts <- log2(filtCounts + 1)
boxplot(logCounts, main = 'Log2 counts', las = 2)

rlogcounts <- rlog(filtCounts)
boxplot(rlogcounts, main = 'rlog counts', las = 2)

pcDat <- prcomp(t(rlogcounts))
autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5)

?autoplot.prcomp

#exercise
autoplot(pcDat,
         data = sampleinfo,
         x = 2,
         y = 3,
         colour = "Status",
         shape = "TimePoint",
         size = 5)

library(ggrepel)

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5) + 
  geom_text_repel(aes(x = PC1, y = PC2, label = SampleName), box.padding = 0.8)


sampleinfo <- mutate(sampleinfo, Status = case_when(
  SampleName=="SRR7657882" ~ "Uninfected",
  SampleName=="SRR7657873" ~ "Infected",
  TRUE ~ Status))

dir.create("results")
write_tsv(sampleinfo, "results/SampleInfo_Corrected.txt")

autoplot(pcDat,
         data = sampleinfo,
         colour = "Status",
         shape = "TimePoint",
         size = 5)
