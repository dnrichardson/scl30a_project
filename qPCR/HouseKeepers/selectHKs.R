
library(ReadqPCR)
library(NormqPCR)
library(dplyr)

qPCRBatch.qPCR <- read.qPCR("qPCR.example.txt")

test <- read.qPCR("test.txt")

tmp <- read.table("test2.txt", header = TRUE) %>% select(Sample:Cq) %>%
        filter(Sample != "NTC")

write.table(tmp, file = "cleaned.txt", quote = FALSE, row.names = FALSE)

test2 <- read.qPCR("cleaned.txt")

rownames(exprs(test))
exprs(test)
exprs.well.order(test)

str(exprs(test))

## Combine the "technical" replicates
combinedData <- combineTechReps(test)
combined2 <- combineTechReps(test2)

## Select most stable housekeepers
resHK <- selectHKs(combined2, method = "geNorm", minNrHKs = 2, log = FALSE, 
                   Symbols = featureNames(combined2))
