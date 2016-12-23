## October 25, 2016
## R script for post-processing output files of rMATS, SplAdder and SUPPA.

## Load in the helper functions
source("processing_functions.R")

library(dplyr)
library(ggplot2)
library(VennDiagram)
library(plyr)
library(tidyr)

## Set basic parameters
dir_prefix <- "~/Dropbox/Current/newscl30a/"
fdr <- 0.05
dpsi <- 0.2
dcut <- 2 ## set to 2 if you want to disable filtering

## Get entire list of gene IDs from ATRT2D in order to query thalemine for gene symbols
## cat AtRTD2_19April2016.gtf | gawk 'BEGIN {FS = ";"} {print $2}' | uniq | gawk '{print $2}' | 
## sed 's/\"//g' > thalemine.genelist.txt

## The following tsv file was obtained from Araport by uploading the gene IDs 
## found in the AtRTD2 April 2016 GTF file. It's gene symbols and descriptions
## will be joined to the output tables from each program
geneNames3 <- tbl_df(read.table("allthalemine.genes.txt.tsv", 
                                stringsAsFactors = FALSE, header = FALSE,
                                sep = "\t", na.strings = "", quote = "", comment.char = "")) %>%
        dplyr::rename(GeneID = V1, symbol = V2, description = V3, is_obsolete = V4)


## Create lists for applying various functions using Sys.glob
rmatsList <- Sys.glob(paste0(dir_prefix, 
                             "atrt2d_MATS/MATS_output/*.MATS.ReadsOnTargetAndJunctionCounts.txt"))
suppaList <- Sys.glob(paste0(dir_prefix, "suppa/atrtd2_quasi/*.dpsi"))
suppaPsiVecList <- Sys.glob(paste0(dir_prefix, "suppa/atrtd2_quasi/*.psivec"))
spladderList_mgraphs <- Sys.glob(paste0(dir_prefix, "spladder/*.confirmed.txt"))
spladderList_testR <- Sys.glob(paste0(dir_prefix, "spladder/*.tsv"))

## Re-analyze suppa with TPM > 5 cutoff
suppaList2 <- Sys.glob(paste0(dir_prefix, "suppa/atrtd2_F5/*.dpsi"))
suppaPsiVecList2 <- Sys.glob(paste0(dir_prefix, "suppa/atrtd2_F5/*.psivec"))
## Run the various processing functions on each list

#####################
## RMATS PROCESSING #
#####################

## Get significant events at the event level and get unique genes
myRmatsEvents <- lapply(rmatsList, processRMATS, fdr, dpsi, dcut = 2)
myRmatsGenes <- lapply(rmatsList, processRMATS, fdr, dpsi, dcut = 2, events = FALSE)

## get table of all events with coordinates
myRmatsCoordsList <- lapply(rmatsList, rMATsCoords, fdr, dpsi, dcut = 2)

## Create output table ranked by deltaPSI
paulaTableEventsR <- bind_rows(myRmatsEvents) %>% dplyr::arrange(desc(abs(IncLevelDifference)))
finalTableEventsR <- left_join(paulaTableEventsR, geneNames3, by = "GeneID")
## Get the same table, but with event coordinates
dcut2EventsR <- bind_rows(myRmatsCoordsList)
write.table(finalTableEventsR, 
            file = paste0("AllEvents_ranked_by_dPSI_rMATS_dcut_", 
                          dcut, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)


#########################
## Spladder processing #
########################
mySpladderEvents <- mapply(processSpladder, spladderList_mgraphs, 
                           spladderList_testR, dpsi = dpsi, dcut = dcut, SIMPLIFY = FALSE)
mySpladderGenes <- mapply(processSpladder, spladderList_mgraphs, 
                          spladderList_testR, dpsi = dpsi, dcut = dcut, events = FALSE, 
                          SIMPLIFY = FALSE)

mySpladderCoordList <- mapply(SpladderCoords, spladderList_mgraphs, 
                              spladderList_testR, dpsi = dpsi, dcut = 2, SIMPLIFY = FALSE)

## Create a table arranged by dPSI

paulaTableEventsSp <- bind_rows(mySpladderEvents) %>% dplyr::arrange(desc(abs(dPSI)))
paulaTableEventsSpAll <- bind_rows(mySpladderCoordList) %>% dplyr::arrange(desc(abs(dPSI)))

finalTableEventsSp <- left_join(paulaTableEventsSp, geneNames3, by = "GeneID")
finalTableEventsSpAll <- left_join(paulaTableEventsSpAll, geneNames3, by = "GeneID")
write.table(finalTableEventsSp, 
            file = paste0("AllEvents_ranked_by_dPSI_spladder_dcut_", 
                          dcut, ".txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE)
write.table(finalTableEventsSpAll, 
            file = paste0("AllEvents_ranked_by_dPSI_spladder_dcut_", 
                          dcut, ".txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE)


#####################
## Suppa processing #
#####################
mySuppaEvents <- mapply(processSuppa, suppaList, suppaPsiVecList, fdr, dpsi, dcut, SIMPLIFY = FALSE)
mySuppaGenes <- mapply(processSuppa, suppaList, suppaPsiVecList, fdr, dpsi, dcut, events = FALSE, 
                       SIMPLIFY = FALSE)

## Process the TPM > 5 files -- only fails on MX.. why? Because there are no
## pvalues below 0.05! Therefore, we will exclude MX from the list.

mySuppaEvents2 <- mapply(processSuppa, suppaList2, suppaPsiVecList2, fdr, dpsi, dcut, SIMPLIFY = FALSE)

## Create a table arranged by dPSI
paulaTableEventsSu <- bind_rows(mySuppaEvents) %>% dplyr::arrange(desc(abs(dPSI)))
finalTableEventsSu <- left_join(paulaTableEventsSu, geneNames3, by = "GeneID")
write.table(finalTableEventsSu, 
            file = paste0("AllEvents_ranked_by_dPSI_SUPPA_", 
                          dcut, ".txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE)
## For the TPM > 5
paulaTableEventsSu2 <- bind_rows(mySuppaEvents2) %>% dplyr::arrange(desc(abs(dPSI)))
finalTableEventsSu2 <- left_join(paulaTableEventsSu2, geneNames3, by = "GeneID")
# only 25 events at this threshold. Maybe need to try TPM > 2

##########################################################################
## Venn diagram of unique genes across programs                          #
## Need to create a dataframe with geneIDs as variables for each program #
##########################################################################
head(paulaTableEvents)

forVennRmats <- finalTableEventsR$GeneID
forVennSuppa <- finalTableEventsSu$GeneID        
forVennSpladder <- finalTableEventsSp$GeneID

## Create the venn diagram 
venn.diagram(list(rMATS=forVennRmats, SUPPA=forVennSuppa, SplAdder=forVennSpladder), 
             filename = "3way_genes_updated.png",
             imagetype = "png", main = "Overlapping genes across programs", 
             sub = "All Event Types, no filtering on psi of replicates",
             col = "transparent",
             fill = c("cornflowerblue", "green", "darkorchid1"), force.unique = TRUE)

commonAcross3 <- Reduce(intersect, list(forVennRmats, forVennSuppa, forVennSpladder))
rmats_suppa <- intersect(forVennRmats, forVennSuppa)
rmats_spladder <- intersect(forVennRmats, forVennSpladder)
suppa_spladder <- intersect(forVennSuppa, forVennSpladder)

## Overlap of rmats and suppa (368) - 10 in common with all 3 programs

onlyrMats_suppa <- setdiff(rmats_suppa, commonAcross3)
onlyrMats_spladder <- setdiff(rmats_spladder, commonAcross3)

## subset AllEvents ranked by dpsi for rMATS of the 358 genes
## Note, you will get more than 358!

onlyrMats_suppa_table <- finalTableEventsR[finalTableEventsR$GeneID %in% onlyrMats_suppa, ]

## get the suppa version of the table
onlySuppa_rmats_table <- finalTableEventsSu[finalTableEventsSu$GeneID %in% onlyrMats_suppa, ]

##onlyrMats_spladder <- finalTableEventsR[finalTableEventsR$GeneID %in% onlyrMats_suppa,]
write.table(onlyrMats_suppa_table, 
            file = paste0("AllEvents_ranked_by_dPSI_only_rmats_suppa_dcut_", 
                          dcut, ".txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE)
write.table(onlySuppa_rmats_table, 
            file = paste0("AllEvents_ranked_by_dPSI_only_suppa_rmats_dcut_", 
                          dcut, ".txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE)

## Sanity check, unique genes should be 358
length(unique(onlyrMats_suppa_table$GeneID))
length(unique(onlySuppa_rmats_table$GeneID))

## With accompanying table with only geneIDs, have to join by hand to make excel table of uneven rows

write.table(commonAcross3, file = paste0("commonAcross3_dcut_", dcut, ".txt"), quote = FALSE, 
            row.names = FALSE)
write.table(rmats_suppa, file = paste0("rmats_suppa_dcut_", dcut, ".txt"), quote = FALSE, 
            row.names = FALSE)
write.table(rmats_spladder, file = paste0("rmats_spladder_dcut_", dcut, ".txt"), quote = FALSE, 
            row.names = FALSE)
write.table(suppa_spladder, file = paste0("suppa_spladder_dcut_", dcut, ".txt"), quote = FALSE, 
            row.names = FALSE)

## Pull out the genes common across all 3 programs, with each program's output

common10.rMATS <- finalTableEventsR[finalTableEventsR$GeneID %in% commonAcross3,]
common10.spladder <- finalTableEventsSp[finalTableEventsSp$GeneID %in% commonAcross3,]
common10.suppa <- finalTableEventsSu[finalTableEventsSu$GeneID %in% commonAcross3,]

write.table(common10.rMATS, 
            file = paste0("AllEvents_ranked_by_dPSI_common10.rmats_dcut_", dcut, ".txt"), 
            quote = FALSE, sep = "\t", 
            row.names = FALSE)

write.table(common10.spladder, 
            file = paste0("AllEvents_ranked_by_dPSI_common10.spladder_dcut_", 
                          dcut, ".txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE)

write.table(common10.suppa, 
            file = paste0("AllEvents_ranked_by_dPSI_common10.suppa_dcut_", 
                          dcut, ".txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE)


#########################################################################
## Create barplots of AS event counts for each program alongside unique #
##  gene counts for that event                                          #
#########################################################################

## Order for SplAdder is: A3, A5, SE, RI, MSE
Types <- c("A3", "A5", "SE", "RI", "MSE")
Counts <- sapply(mySpladderEvents, nrow)
uniqGeneCounts <- sapply(mySpladderGenes, nrow) 

eventsDf <- tbl_df(data.frame(Types, Counts)) %>% 
        dplyr::rename(types = Types, counts = Counts)
eventsDf$variable <- "events"

genesDf <- tbl_df(data.frame(Types, uniqGeneCounts)) %>% 
        dplyr::rename(types = Types, counts = uniqGeneCounts)
genesDf$variable <- "genes"

combinedDf <- bind_rows(eventsDf, genesDf)

## Create plot of AS event types

splad_ASCounts <- ggplot(combinedDf, aes(x = types, y = counts, fill = variable)) + 
        geom_bar(stat = "identity", position = "dodge") + 
        ggtitle(paste0("splAdder counts of AS event types (dPSI > 0.2)\n total counts of AS events: ",
sum(Counts), ", AS genes: ", sum(uniqGeneCounts))) +
        labs( x = "AS event types" ) + geom_text(aes(label=counts), 
                                                 position = position_dodge(width = 1), vjust = 0) +
        guides(fill=guide_legend(title=NULL)) + theme_grey()

ggsave(paste0("SplAddersASCounts_dcut_", dcut, ".png"), splad_ASCounts)

## rMATS
## Order is: A3, A5, MXE, RI, SE
Types <- c("A3", "A5", "MXE", "RI", "SE")
Counts <- sapply(myRmatsEvents, nrow)
uniqGeneCounts <- sapply(myRmatsGenes, nrow)

eventsDf <- tbl_df(data.frame(Types, Counts)) %>% 
        dplyr::rename(types = Types, counts = Counts)
eventsDf$variable <- "events"

genesDf <- tbl_df(data.frame(Types, uniqGeneCounts)) %>% 
        dplyr::rename(types = Types, counts = uniqGeneCounts)
genesDf$variable <- "genes"

combinedDf <- bind_rows(eventsDf, genesDf)

## Create plot of AS event types

rmatsASCounts <- ggplot(combinedDf, aes(x = types, y = counts, fill = variable)) + 
        geom_bar(stat = "identity", position = "dodge") + 
        ggtitle(paste0("rMATS counts of AS event types (dPSI > 0.2)\n total counts of AS events: ",
                       sum(Counts), ", AS genes: ", sum(uniqGeneCounts))) +
        labs( x = "AS event types" ) + geom_text(aes(label=counts), 
                                                 position = position_dodge(width = 1), vjust = 0) +
        guides(fill=guide_legend(title=NULL)) + theme_grey()

ggsave(paste0("rMATSsASCounts_dcut_", dcut, ".png"), rmatsASCounts)

## SUPPA
## Order is: A3, A5, MXE, RI, SE
Types <- c("A3", "A5", "MXE", "RI", "SE")
Counts <- sapply(mySuppaEvents, nrow)
uniqGeneCounts <- sapply(mySuppaGenes, nrow)

eventsDf <- tbl_df(data.frame(Types, Counts)) %>% 
        dplyr::rename(types = Types, counts = Counts)
eventsDf$variable <- "events"

genesDf <- tbl_df(data.frame(Types, uniqGeneCounts)) %>% 
        dplyr::rename(types = Types, counts = uniqGeneCounts)
genesDf$variable <- "genes"

combinedDf <- bind_rows(eventsDf, genesDf)

## Create plot of AS event types

SuppaASCounts <- ggplot(combinedDf, aes(x = types, y = counts, fill = variable)) + 
        geom_bar(stat = "identity", position = "dodge") + 
        ggtitle(paste0("SUPPA counts of AS event types (dPSI > 0.2)\n total counts of AS events: ",
                       sum(Counts), ", AS genes: ", sum(uniqGeneCounts))) +
        labs( x = "AS event types" ) + 
        geom_text(aes(label=counts), position = position_dodge(width = 1), vjust = 0) +
        guides(fill=guide_legend(title=NULL)) + theme_grey()

ggsave(paste0("SUPPAsASCounts_dcut_", dcut, ".png"), SuppaASCounts)


## not happy with the below
grid.arrange(splad_ASCounts, rmatsASCounts, SuppaASCounts, ncol = 3, heights = 1:2)

## Look at distributions of differences between replicates.
## Looks like a good cutoff would be 0.3

png("histogram_diff_psi_reps.png")

par(mfrow = c(3,2))

h <- hist(finalTableEventsR$dWT, plot = FALSE)
h$density <- h$counts/sum(h$counts)*100
plot(h, freq = FALSE, ylab = "percent", xlab = "difference in PSI between reps", 
     main = "Distribution of differences\n in WT PSI for rMATS events")
abline(v = 0.2)
#rug(h$counts)

h <- hist(finalTableEventsR$dKO, plot = FALSE)
h$density <- h$counts/sum(h$counts)*100
plot(h, freq = FALSE, ylab = "percent", xlab = "difference in PSI between reps", 
     main = "Distribution of differences\n in KO PSI for rMATS events")
abline(v = 0.2)


h <- hist(finalTableEventsSp$dWT, plot = FALSE)
h$density <- h$counts/sum(h$counts)*100
plot(h, freq = FALSE, ylab = "percent", xlab = "difference in PSI between reps", 
     main = "Distribution of differences\n in WT PSI for SplAdder events")
abline(v = 0.2)

h <- hist(finalTableEventsSp$dKO, plot = FALSE)
h$density <- h$counts/sum(h$counts)*100
plot(h, freq = FALSE, ylab = "percent", xlab = "difference in PSI between reps", 
     main = "Distribution of differences\n in KO PSI for SplAdder events")
abline(v = 0.2)


h <- hist(finalTableEventsSu$dWT, plot = FALSE)
h$density <- h$counts/sum(h$counts)*100
plot(h, freq = FALSE, ylab = "percent", xlab = "difference in PSI between reps", 
     main = "Distribution of differences\n in WT PSI for SUPPA events")
abline(v = 0.2)


h <- hist(finalTableEventsSu$dKO, plot = FALSE)
h$density <- h$counts/sum(h$counts)*100
plot(h, freq = FALSE, ylab = "percent", xlab = "difference in PSI between reps", 
     main = "Distribution of differences\n in KO PSI for SUPPA events")
abline(v = 0.2)

dev.off()

par(mfrow=c(1,1))

## Plot the PSI values against each other
## rMATS
plot(finalTableEventsR$atRWT2S.psi, finalTableEventsR$atRWT3S.psi, 
     col = rgb(0, 0, 0, 0.2), pch = 19)
abline(lm(finalTableEventsR$atRWT3S.psi~finalTableEventsR$atRWT2S.psi), 
       col="red") # regression line (y~x) 
cor(finalTableEventsR$atRWT2S.psi, finalTableEventsR$atRWT3S.psi)

plot(finalTableEventsR$atRKO1S.psi, finalTableEventsR$atRKO3S.psi, 
     col = rgb(0, 0, 0, 0.2), pch = 19)

## SplAdder
plot(finalTableEventsSp$atRWT2S.psi, finalTableEventsSp$atRWT3S.psi, 
     col = rgb(0, 0, 0, 0.2), pch = 19)
plot(finalTableEventsSp$atRKO1S.psi, finalTableEventsSp$atRKO3S.psi, 
     col = rgb(0, 0, 0, 0.2), pch = 19)
cor(finalTableEventsSp$atRWT2S.psi, finalTableEventsSp$atRWT3S.psi)
cor(finalTableEventsSp$atRKO1S.psi, finalTableEventsSp$atRKO3S.psi)

## SUPPA
plot(finalTableEventsSu$atRWT2S.psi, finalTableEventsSu$atRWT3S.psi, 
     col = rgb(0, 0, 0, 0.2), pch = 19)
plot(finalTableEventsSu$atRKO1S.psi, finalTableEventsSu$atRKO3S.psi, 
     col = rgb(0, 0, 0, 0.2), pch = 19)
cor(finalTableEventsSu$atRWT2S.psi, finalTableEventsSu$atRWT3S.psi)
cor(finalTableEventsSu$atRKO1S.psi, finalTableEventsSu$atRKO3S.psi)

## Convert Type column to factor
finalTableEventsR <- transform(finalTableEventsR, Type = factor(Type))
table(finalTableEventsR$Type)

finalTableEventsSu <- transform(finalTableEventsSu, Type = factor(Type))
table(finalTableEventsSu$Type)

finalTableEventsSp <- transform(finalTableEventsSp, Type = factor(Type))
table(finalTableEventsSp$Type)

## Use GRanges to see which events actually overlap each other

## Create a dataframe that holds the coordinates of the AS event in question. Lets try SplAdder IR for now
## merge_graphs_intron_retention_C3.confirmed.txt
## contig	strand	event_id	gene_name	exon1_start	exon1_end	intron_start	intron_end	exon2_start	exon2_end
## Chr2	+	intron_retention_42335	AT2G43010 17887068	17887157

##########################
#       SPLADDER         #
##########################

## Build the SplAdder data frame from the output of the function, SpladderCoords

### RI ###
spladderRI.df <- mySpladderCoordList[[4]] %>% 
        dplyr::select(contig:gene_name, intron_start, intron_end) %>% 
        dplyr::rename(
        chr = contig, start = intron_start, end = intron_end, GeneID = gene_name)

grSpladderRI <- makeGRangesFromDataFrame(spladderRI.df, keep.extra.columns = TRUE)

### A3 ###
spladderA3_1.df <- mySpladderCoordList[[1]] %>% 
        dplyr::select(contig:gene_name, exon_alt1_start, exon_alt1_end) %>% 
        dplyr::rename(
        chr = contig, start = exon_alt1_start, end = exon_alt1_end, GeneID = gene_name)

spladderA3_2.df <- mySpladderCoordList[[1]] %>% 
        dplyr::select(contig:gene_name, exon_alt2_start, exon_alt2_end) %>% 
        dplyr::rename(
        chr = contig, start = exon_alt2_start, end = exon_alt2_end, GeneID = gene_name)

grSpladderA31 <- makeGRangesFromDataFrame(spladderA3_1.df, keep.extra.columns = TRUE)
grSpladderA32 <- makeGRangesFromDataFrame(spladderA3_2.df, keep.extra.columns = TRUE)

### A5 ###

spladderA5_1.df <- mySpladderCoordList[[2]] %>% 
        dplyr::select(contig:gene_name, exon_alt1_start, exon_alt1_end) %>% 
        dplyr::rename(
        chr = contig, start = exon_alt1_start, end = exon_alt1_end, GeneID = gene_name)

spladderA5_2.df <- mySpladderCoordList[[2]] %>% 
        dplyr::select(contig:gene_name, exon_alt2_start, exon_alt2_end) %>% 
        dplyr::rename(
        chr = contig, start = exon_alt2_start, end = exon_alt2_end, GeneID = gene_name)

grSpladderA51 <- makeGRangesFromDataFrame(spladderA5_1.df, keep.extra.columns = TRUE)
grSpladderA52 <- makeGRangesFromDataFrame(spladderA5_2.df, keep.extra.columns = TRUE)

### ES ###
spladderES.df <- mySpladderCoordList[[3]] %>% 
        dplyr::select(contig:gene_name, exon_start, exon_end) %>% 
        dplyr::rename(
        chr = contig, start = exon_start, end = exon_end, GeneID = gene_name)

grSpladderES <- makeGRangesFromDataFrame(spladderES.df, keep.extra.columns = TRUE)

#########################
#       SUPPA           #
#########################

## Build the gr object from the Suppa Table. "SuAll" means no dcut filter
SuppaRI.df <- finalTableEventsSuAll[finalTableEventsSuAll$Type == "RI",] %>% 
        dplyr::select(event, GeneID)
SuppaA3.df <- finalTableEventsSuAll[finalTableEventsSuAll$Type == "A3",] %>% 
        dplyr::select(event, GeneID)
SuppaA5.df <- finalTableEventsSuAll[finalTableEventsSuAll$Type == "A5",] %>% 
        dplyr::select(event, GeneID)
SuppaES.df <- finalTableEventsSuAll[finalTableEventsSuAll$Type == "SE",] %>% 
        dplyr::select(event, GeneID)
SuppaMXE.df <- finalTableEventsSuAll[finalTableEventsSuAll$Type == "MXE",] %>% 
        dplyr::select(event, GeneID)
nrow(SuppaMXE.df)
# 0 MXE events

## Chromosome
SuppaRI.df$chr <- sapply(strsplit(SuppaRI.df$event, ":"), "[", 2)
SuppaA3.df$chr <- sapply(strsplit(SuppaA3.df$event, ":"), "[", 2)
SuppaA5.df$chr <- sapply(strsplit(SuppaA5.df$event, ":"), "[", 2)
SuppaES.df$chr <- sapply(strsplit(SuppaES.df$event, ":"), "[", 2)

## RI coordinates
tmp <- sapply(strsplit(SuppaRI.df$event, ":"), "[", 4)

SuppaRI.df$start <- sapply(strsplit(tmp, "-"), "[", 1)
SuppaRI.df$end <- sapply(strsplit(tmp, "-"), "[", 2)

SuppaRI.df$strand <- sapply(strsplit(SuppaRI.df$event, ":"), "[", 6)
head(SuppaRI.df)

grSuppaRI <- makeGRangesFromDataFrame(SuppaRI.df, keep.extra.columns = TRUE)

## A3 coordinates ##
a3_1 <- sapply(strsplit(SuppaA3.df$event, ":"), "[", 3)
SuppaA3.df$start <- sapply(strsplit(a3_1, "-"), "[", 1)
SuppaA3.df$end <- sapply(strsplit(a3_1, "-"), "[", 2)
SuppaA3.df$strand <- sapply(strsplit(SuppaA3.df$event, ":"), "[", 5)

## Copy the first A3 dataframe that held info for the first A3 event
## so we can have the coordinates of the second A3 event
SuppaA3.df2 <- dplyr::select(SuppaA3.df, event, GeneID, chr, strand)

a3_2 <- sapply(strsplit(SuppaA3.df$event, ":"), "[", 4)
SuppaA3.df2$start <- sapply(strsplit(a3_1, "-"), "[", 1)
SuppaA3.df2$end <- sapply(strsplit(a3_1, "-"), "[", 2)

grSuppaA3_1 <- makeGRangesFromDataFrame(SuppaA3.df, keep.extra.columns = TRUE)
grSuppaA3_2 <- makeGRangesFromDataFrame(SuppaA3.df2, keep.extra.columns = TRUE)

## A5 coordinates ##

a5_1 <- sapply(strsplit(SuppaA5.df$event, ":"), "[", 3)
SuppaA5.df$start <- sapply(strsplit(a5_1, "-"), "[", 1)
SuppaA5.df$end <- sapply(strsplit(a5_1, "-"), "[", 2)
SuppaA5.df$strand <- sapply(strsplit(SuppaA5.df$event, ":"), "[", 5)

## Copy the first A5 dataframe that held info for the first A5 event 
## so we can also have the coordinates of the second A5 event
SuppaA5.df2 <- dplyr::select(SuppaA5.df, event, GeneID, chr, strand)

a5_2 <- sapply(strsplit(SuppaA5.df$event, ":"), "[", 4)
SuppaA5.df2$start <- sapply(strsplit(a5_1, "-"), "[", 1)
SuppaA5.df2$end <- sapply(strsplit(a5_1, "-"), "[", 2)

grSuppaA5_1 <- makeGRangesFromDataFrame(SuppaA5.df, keep.extra.columns = TRUE)
grSuppaA5_2 <- makeGRangesFromDataFrame(SuppaA5.df2, keep.extra.columns = TRUE)

## SE coordinates

SE <- sapply(strsplit(SuppaES.df$event, ":"), "[", 4)

SuppaES.df$start <- sapply(strsplit(SE, "-"), "[", 1)
SuppaES.df$end <- sapply(strsplit(SE, "-"), "[", 2)

SuppaES.df$strand <- sapply(strsplit(SuppaES.df$event, ":"), "[", 5)
head(SuppaES.df)

grSuppaES <- makeGRangesFromDataFrame(SuppaES.df, keep.extra.columns = TRUE)


###########################
#       RMATS             #
###########################

## RI
rMATSRI.df <- myRmatsCoordsList[[4]] %>% 
        dplyr::select(chr, GeneID, ID, strand, riExonStart_0base, riExonEnd) %>%
        dplyr::rename(start = riExonStart_0base, end = riExonEnd)

nrow(rMATSRI.df)
# 508
grRmatsRI <- makeGRangesFromDataFrame(rMATSRI.df, keep.extra.columns = TRUE)


## A3
rMATSA3_1.df <- myRmatsCoordsList[[1]] %>% 
        dplyr::select(chr, GeneID, ID, strand, longExonStart_0base, longExonEnd) %>%
        dplyr::rename(start = longExonStart_0base, end = longExonEnd)

rMATSA3_2.df <- myRmatsCoordsList[[1]] %>% 
        dplyr::select(chr, GeneID, ID, strand, shortES, shortEE) %>%
        dplyr::rename(start = shortES, end = shortEE)

grRmatsA3_1 <- makeGRangesFromDataFrame(rMATSA3_1.df, keep.extra.columns = TRUE)
grRmatsA3_2 <- makeGRangesFromDataFrame(rMATSA3_2.df, keep.extra.columns = TRUE)

## A5

rMATSA5_1.df <- myRmatsCoordsList[[2]] %>% 
        dplyr::select(chr, GeneID, ID, strand, longExonStart_0base, longExonEnd) %>%
        dplyr::rename(start = longExonStart_0base, end = longExonEnd)

rMATSA5_2.df <- myRmatsCoordsList[[2]] %>% 
        dplyr::select(chr, GeneID, ID, strand, shortES, shortEE) %>%
        dplyr::rename(start = shortES, end = shortEE)


grRmatsA5_1 <- makeGRangesFromDataFrame(rMATSA5_1.df, keep.extra.columns = TRUE)

grRmatsA5_2 <- makeGRangesFromDataFrame(rMATSA5_2.df, keep.extra.columns = TRUE)

## SE ##

rMATSES.df <- myRmatsCoordsList[[5]] %>% 
        dplyr::select(chr, GeneID, ID, strand, exonStart_0base, exonEnd) %>%
        dplyr::rename(start = exonStart_0base, end = exonEnd)

nrow(rMATSES.df)
# 430
grRmatsES <- makeGRangesFromDataFrame(rMATSES.df, keep.extra.columns = TRUE)

#########################################################
# OVERLAP DETECTION OF EVENTS WITH GRANGES ALL PROGRAMS #
#########################################################

### RI ###

## Try with mergeByOverlaps
RI_spl_sup <- mergeByOverlaps(grSpladderRI, grSuppaRI)
# 8 overlaps
RI_spl_rmats <-  mergeByOverlaps(grSpladderRI, grRmatsRI)
# 0
RI_sup_rmats <- mergeByOverlaps(grSuppaRI, grRmatsRI, minoverlap = 10, type = c("any"))
# 197 overlaps -- 185 unique events in Suppa, 169 in rMATS

### A3 ###

## Try with mergeByOverlaps
A3_spl_sup <- mergeByOverlaps(grSpladderA31, grSuppaA3_2, minoverlap = 10)
# 0 overlaps
A3_spl_rmats <-  mergeByOverlaps(grSpladderA31, grRmatsA3_1, minoverlap = 10)
# 5 overlaps
A3_sup_rmats <- mergeByOverlaps(grSuppaA3_1, grRmatsA3_1, minoverlap = 10, type = c("any"))
# 12 overlaps

### A5 ###

## Try with mergeByOverlaps
A5_spl_sup <- mergeByOverlaps(grSpladderA51, grSuppaA5_2, minoverlap = 10)
# 0 overlaps
A5_spl_rmats <-  mergeByOverlaps(grSpladderA51, grRmatsA5_1, minoverlap = 10)
# 4 overlaps
A5_sup_rmats <- mergeByOverlaps(grSuppaA5_1, grRmatsA5_1, minoverlap = 10, type = c("any"))
# 5 overlaps

### ES ###

## Try with mergeByOverlaps
ES_spl_sup <- mergeByOverlaps(grSpladderES, grSuppaES, minoverlap = 10)
# 3 overlaps -- 3 unique events
ES_spl_rmats <-  mergeByOverlaps(grSpladderES, grRmatsES, minoverlap = 10)
# 41 overlaps -- 14 unique events in Spladder, 26 in rMATs
ES_sup_rmats <- mergeByOverlaps(grSuppaES, grRmatsES, minoverlap = 10, type = c("any"))
# 9 overlaps -- 6 unique events in Spladder, 9 in rMATS


#################
# Gviz plotting #
#################

## what I wanted to do was to plot my primer locations on the exons, but
## it is proving to be a real pain in the ass. There is some difficulty in 
## directly plotting my GRanges event object for some reason. To go through
## by hand to plot my primers is a huge time investment that I can't make right
## now. 

options(ucscChromosomeNames=FALSE)

gTrack <- GenomeAxisTrack()


test_event1 <- subset(grSpladderRI, event_id == "intron_retention_42333")

aTrack <- AnnotationTrack(test_event1, shape = "box")

testtx <- GeneRegionTrack(txdb, chromosome = 2, start = 17886000, end = 17890000,
                          transcriptAnnotation = TRUE, min.height = 5)
testtx2 <- GeneRegionTrack(txdb.arab, chromosome = "Chr2", start = 17885000, end = 17890500,
                           transcriptAnnotation = TRUE)

annotTest <- AnnotationTrack(start = 17886651, end = 17886993,
                             chromosome = 2, strand = "+",
                             id = "intron_retention_42333",
                             name = "event", 
                             shape = "box",
                             featureAnnotation = "id",
                             cex.feature = 0.6,
                             fontcolor.feature = 1)


plotTracks(list(gTrack, testtx2))
plotTracks(list(gTrack, testtx, annotTest), extend.left = 500, extend.right = 500)

plotTracks(list(gTrack, test3))

## Try with ggbio
library(ggbio)
library(biovizBase)

xscripts <- transcripts(txdb)
xscripts2 <- transcripts(txdb.arab)

exons2 <- exonsBy(txdb.arab, by = "tx")
wh <- subset(xscripts, tx_name == "AT2G43010.1")
whE <- subsetByOverlaps(exons2, wh2)

tmp <- subset(xscripts2, tx_name = "AT5G11870_c1")
subsetByOverlaps(exons2, tmp)

## Spladder RI Event 42333
test_event1 <- subset(grSpladderRI, event_id == "intron_retention_42333")
wh2 <- subset(xscripts2, tx_name == "AT2G43010_c2")
autoplot(txdb.arab, which = wh2)

# alternatively

one <- ggplot() + geom_alignment(txdb.arab, which = wh2) + 
        geom_rect(data = test_event1, aes(xmin = start(test_event1), xmax = end(test_event1), 
                                          ymin = -Inf, ymax = Inf), 
                      alpha = 0.4, fill = "green"  )

two <- ggplot() + geom_rect(data = test_event1, aes(xmin = start(test_event1), 
                                                    xmax = end(test_event1), 
                                             ymin = -Inf, ymax = Inf), 
                     alpha = 0.4, fill = "green"  )

tracks(one, two)


## Load BSgenome Athaliana
library(BSgenome)
source("http://bioconductor.org/biocLite.R")
biocLite("BSgenome.Athaliana.TAIR.04232008")
library(BSgenome.Athaliana.TAIR.04232008)
genome <- BSgenome.Athaliana.TAIR.04232008

mir162a <- getSeq(genome, "chr5", start = 2633919, end = 2635578, strand = "-" )
as.character(mir162a)

## Get tranposase sequence for PCR primer design, which includes the 
## exonic sequences before and after the retained intron
transposase <- mySpladderCoordList[[4]] %>% 
        filter(event_id == "intron_retention_102273") %>% 
        dplyr::select(contig:exon2_end)

transposase.df <- dplyr::rename(transposase, chr = contig)

mydf <- data.frame(chr = rep("Chr5", 3), start = c(12664506, 12664660, 12664809),
                   end = c(12664659, 12664808, 12665035), 
                   event_id = rep("intron_retention_102273"),
                   strand = rep("+",3))

grmydf <- makeGRangesFromDataFrame(mydf, keep.extra.columns = TRUE)
transposase.exonPre_exonPost <- getSeq(genome, "chr5", 
                                       start = 12664506, end = 12665035,
                                       strand = "+")

tranposase.exonPre <- getSeq(genome, "chr5", start = 12664506, 
                             end = 12664659, strand = "+")
transposase.exonPost <- getSeq(genome, "chr5", start = 12664809, 
                               end = 12665035, strand = "+")

paste0(as.character(tranposase.exonPre), as.character(transposase.exonPost))

transposase.genomic <- as.character(transposase.exonPre_exonPost)

## Read in list of DEGs and check for how many overlap with the unique
## AS GeneIDs for each program

DEGs <- read.csv("../QoRTS_AtRT2D/KO_vs_WT.csv", stringsAsFactors = FALSE) 

intersect(DEGs$ID, forVennSpladder)
# [1] "AT5G33395" "AT3G30720" "AT4G28520"
DEGs[DEGs$ID %in% intersect(DEGs$ID, forVennSpladder),]
# ID    baseMean log2FoldChange     lfcSE      stat                          pvalue                        padj    Symbol               Description
# 3  AT5G33395    145.1317      -1.924632 0.2704467 -7.116492 0.00000000000110999999999999992 0.0000000004030000000000000 AT5G33395 transposable_element_gene
# 5  AT3G30720    161.4574      -1.821319 0.2855626 -6.378003 0.00000000017900000000000000226 0.0000000333000000000000004       QQS          qua-quine starch
# 44 AT4G28520 242071.9252       1.933505 0.1801575 10.732299 0.00000000000000000000000000718 0.0000000000000000000000636      CRU3              cruciferin 3
intersect(DEGs$ID, forVennSuppa)
# [1] "AT1G43590" "AT1G51680"
DEGs[DEGs$ID %in% intersect(DEGs$ID, forVennSuppa),]
# ID  baseMean log2FoldChange     lfcSE     stat                   pvalue                   padj    Symbol               Description
# 26 AT1G43590  136.5881       1.610701 0.2745619 5.866439 0.0000000044500000000000 0.00000051900000000000 AT1G43590 transposable_element_gene
# 48 AT1G51680 1276.7257       2.006118 0.2074360 9.671023 0.0000000000000000000004 0.00000000000000000118      4CL1  4-coumarate:CoA ligase 1
intersect(DEGs$ID, forVennRmats)
# [1] "AT3G13570" "AT1G43590" "AT3G15450"
DEGs[DEGs$ID %in% intersect(DEGs$ID, forVennRmats),]
# ID   baseMean log2FoldChange     lfcSE      stat                     pvalue                    padj    Symbol
# 2  AT3G13570  1333.4100      -2.057119 0.2101707 -9.787847 0.000000000000000000000127 0.000000000000000000553    SCL30A
# 26 AT1G43590   136.5881       1.610701 0.2745619  5.866439 0.000000004450000000000000 0.000000519000000000000 AT1G43590
# 43 AT3G15450 20298.5906       1.862172 0.2073780  8.979603 0.000000000000000000272000 0.000000000000000425000 AT3G15450
# Description
# 2                      SC35-like splicing factor 30A
# 26                         transposable_element_gene
# 3 aluminum induced protein with YGL and LRDR motifs


## Some test filtering of the rMATS output to only have a certain number of Inclusion level counts
## so as to remove potential false positives

thresh <- 50

separatedRmats <- dcut2EventsR %>% 
        tidyr::separate(IC_SAMPLE_1, into = c("wtIC1", "wtIC2"), sep = ",") %>%
        tidyr::separate(IC_SAMPLE_2, into = c("koIC1", "koIC2"), sep = ",") %>%
        tidyr::separate(SC_SAMPLE_1, into = c("wtSC1", "wtSC2"), sep = ",") %>%
        tidyr::separate(SC_SAMPLE_2, into = c("koSC1", "koSC2"), sep = ",")

separatedRmats %>% filter(wtIC1 > thresh, wtIC2 > thresh, koIC1 > thresh, koIC2 > thresh,
                          wtSC1 > thresh, wtSC2 > thresh, koSC1 > thresh, koSC2 > thresh)

## Higher inclusion levels in WT, higher skipping levels in KO
separatedRmats %>% filter(wtIC1 > thresh, wtIC2 > thresh, koSC1 > thresh, koSC2 > thresh) %>%
        select(ID, GeneID, Type)

## Higher skipping levels in WT, higher inclusion levels in KO
separatedRmats %>% filter(wtSC1 > thresh, wtSC2 > thresh, koIC1 > thresh, koIC2 > thresh) %>%
        select(ID, GeneID, Type)

## Load up the DeSeq2 DDS Normalized Gene counts "allNormalizedCounts"

load("deseq2NormCounts.Rda")

str(allNormalizedCounts)
class(allNormalizedCounts)

## It's a matrix. Need to convert to DF. 
## Remove scientific notation, to add back, scipen = 0
options(scipen=999)
allNormalizedCounts <- as.data.frame(allNormalizedCounts)

## Use these normalized expression counts to filter out "lowly expressed genes" in my
## AS outputs

## Convert row names to variable, 
library(tibble)
allNormalizedCounts <- rownames_to_column(allNormalizedCounts)
normCounts <- allNormalizedCounts %>% tbl_df() %>% mutate(GeneID = rowname) %>%
        select(GeneID, atRWT2S:atRKO3S)
normCounts$baseMean <- rowMeans(normCounts[c(2:5)])

## Log transform baseMean and plot histogram

normCounts$logBaseMean <- log10(normCounts$baseMean)

ggplot(normCounts, aes(x=logBaseMean)) +
        geom_histogram(binwidth=.25, colour="black", fill="white") +
        geom_vline(xintercept = 2, colour = "green") +
        labs(title = "Distribution of DeSeq2 log transformed normalized counts")

## bimodal distribution 

## try plotly
library(plotly)

plot_ly(alpha = 0.3) %>%
        add_histogram(x = ~log(normCounts$atRWT2S)) %>%
        add_histogram(x = ~log(normCounts$atRWT3S)) %>%
        add_histogram(x = ~log(normCounts$atRKO1S)) %>%
        add_histogram(x = ~log(normCounts$atRKO3S)) %>%
        layout(barmode = "overlay")

## Join to each AS's program finalTable

finalSpladderGeneExpCounts <- left_join(finalTableEventsSpAll, normCounts, by = "GeneID")
finalSuppaGeneExpCounts <- left_join(finalTableEventsSu, normCounts, by = "GeneID")
finalRMatsGeneExpCounts <- left_join(separatedRmats, normCounts, by = "GeneID")

## Plot logBaseMean vs dPSI for spladder

tmp <- data.frame(x = finalSpladderGeneExpCounts$logBaseMean, y = abs(finalSpladderGeneExpCounts$dPSI))

ggplot(tmp, aes(x = x, y = y ) ) + geom_point() +
        labs(title = "logBaseMean versus dPSI for Spladder Events", 
             x = "Log Base Mean", y = "absolulte dPSI") + geom_smooth()

## now for suppa
tmp <- data.frame(x = finalSuppaGeneExpCounts$logBaseMean, y = abs(finalSuppaGeneExpCounts$dPSI))

ggplot(tmp, aes(x = x, y = y ) ) + geom_point() +
        labs(title = "logBaseMean versus dPSI for Suppa Events", 
             x = "Log Base Mean", y = "absolulte dPSI") + geom_smooth()

## now for rmats
tmp <- data.frame(x = finalRMatsGeneExpCounts$logBaseMean, y = abs(finalRMatsGeneExpCounts$IncLevelDifference))
#tmp2 <- data.frame(x = finalRMatsGeneExpCounts$baseMean, y = abs(finalRMatsGeneExpCounts$IncLevelDifference))
ggplot(tmp, aes(x = x, y = y ) ) + geom_point() +
        labs(title = "logBaseMean versus dPSI for rMATS Events", 
             x = "Log Base Mean", y = "absolulte dPSI") + geom_smooth()

#ggplot(tmp2, aes(x = x, y = y ) ) + geom_point() +
#       labs(title = "BaseMean versus dPSI for rMATS Events", 
#             x = "Base Mean", y = "absolulte dPSI") + geom_smooth()

## Now, filter the tables based on gene expression. The gene should be moderately 
## expressed in at least one condition! Perhaps better to filter such that each replicate should
## be expressed at greater than threshold.

basemean.thresh <- 500
## log of basemean threshold
logbm <- 2

filteredSpladder <- finalSpladderGeneExpCounts %>% 
        filter(baseMean > basemean.thresh) %>% 
        filter(baseMean < 1000)

## Get retained introns
filteredSpladderRI <- finalSpladderGeneExpCounts %>% 
        filter(logBaseMean > logbm) %>%
        filter(grepl("intron_retention", event_id))

## Get skipped exons
filteredSpladderSE <- finalSpladderGeneExpCounts %>% 
        filter(logBaseMean > logbm) %>%
        filter(grepl("exon_skip", event_id))


## remove NA columns
filteredSpladderRI <- filteredSpladderRI[colSums(!is.na(filteredSpladderRI)) > 0]
filteredSpladderSE <- filteredSpladderSE[colSums(!is.na(filteredSpladderSE)) > 0]

## consider only genes that have an intron coverage difference > 
## intron_cov: mean coverage of the retained intron

filteredSpladderRI$absCovDiff <- 
        abs((filteredSpladderRI$atRWT2S.intron_cov + filteredSpladderRI$atRWT3S.intron_cov) / 2 -
                    (filteredSpladderRI$atRKO1S.intron_cov + filteredSpladderRI$atRKO3S.intron_cov) / 2)

SpladderRI.20 <- filteredSpladderRI %>% filter(absCovDiff > 20)
SpladderRI.15 <- filteredSpladderRI %>% filter(absCovDiff > 15)
SpladderRI.10 <- filteredSpladderRI %>% filter(absCovDiff > 10)

plot(filteredSpladderRI$logBaseMean, filteredSpladderRI$dPSI)
hist(filteredSpladder$baseMean, bins )

filteredSuppa <- finalSuppaGeneExpCounts %>%
        filter(baseMean > basemean.thresh) 

filteredRmats <- finalRMatsGeneExpCounts %>%
        filter(baseMean > basemean.thresh) %>%
        arrange(desc(baseMean, IncLevelDifference)) %>%
        select(ID, GeneID, Type, IncLevelDifference, baseMean) %>%
        filter(Type == "RI")

filteredRmatsRI <- finalRMatsGeneExpCounts %>%
        filter(logBaseMean > 2.5) %>%
        arrange(desc(IncLevelDifference)) %>%
        select(ID, GeneID, Type, IncLevelDifference, baseMean) %>%
        filter(Type == "RI")


## The above is nice, but I think I need to somehow involve junction read counts

## jSplice results ############################################################
## Read in jSplice data that has been pre-processed by hand to remove spurious, 
## aggregated ASMs
## cat tmp.50.txt | grep -v "Unknown"| awk 'BEGIN {FS = "|"}; {print $1"\t"$2}' 
##| cut -f 1,3,4,5,6,7,8 > jSplice.cleaned.inc50.txt
## Load this file into excel, export as .csv

## Read in the inclusion level = 0.5 set
jSplice <- read.csv("../tmp/jSplice.cleaned.inc50.csv", 
                      stringsAsFactors = FALSE)

## Filter for genes that have fdr > 0.10 in both replicates
jSpliceFiltered <- jSplice %>% filter(fdrRepA < 0.10 & fdrRepB < 0.10)

## Create vector of geneIDs and see overlap with other programs
jSpliceVec <- jSpliceFiltered$GeneID

sum(jSpliceVec %in% finalSpladderGeneExpCounts$GeneID)
# 2 in common with spladder
jSpliceVec[jSpliceVec %in% finalSpladderGeneExpCounts$GeneID]
# [1] "AT3G02300" "AT1G79245"

## rMATS
sum(jSpliceVec %in% finalRMatsGeneExpCounts$GeneID)
# 8 in common with rMATs, 2 of which are the same as above
jSpliceVec[jSpliceVec %in% finalRMatsGeneExpCounts$GeneID]
# [1] "AT3G59430" "AT2G17730" "AT3G02300" "AT5G65685" "AT1G18750" "AT1G79245" "AT3G06210"
# [8] "AT1G27630"

## SUPPA
sum(jSpliceVec %in% finalSuppaGeneExpCounts$GeneID)
# 11
jSpliceVec[jSpliceVec %in% finalSuppaGeneExpCounts$GeneID]
# [1] "AT3G59430" "AT4G01550" "AT4G30570" "AT1G18750" "AT1G79245" "AT1G63670" "AT3G06210"
# [8] "AT5G45710" "AT2G28940" "AT1G27630" "AT2G01730"

## Read in 25% inclusion jSplice run
## Filter file on command line:
## cat jSplice_inc25_html_to_txt.txt | awk 'NF > 2' | grep -v "Unknown" | awk 'BEGIN {FS = "|"}; 
## {print $1"\t"$2}' | cut -f 1,3,4,5,6,7,8 > jSplice_inc25.clean.txt
## add variables in excel: 
## GeneID	Type	meanLog2wt.ko	maxLog2wt.ko.B	maxLog2wt.ko.A	FDRB	FDRA

## read in file and filter on FDR 
loc <- "~/Dropbox/Current/newscl30a/scl30a_project/jSplice_Results/jSplice25.csv"
jSplice25 <- tbl_df(read.csv(loc, stringsAsFactors = FALSE)) %>% 
        filter(FDRB < 0.05 & FDRA < 0.05)

jSplice25 <- left_join(jSplice25, geneNames3, by = "GeneID")

## Attach basemean counts
jSplice25 <- left_join(jSplice25, normCounts, by = "GeneID")

## output table for excel
write.table(jSplice25, file = "jSplice25.tsv", quote = FALSE, row.names = FALSE,
            sep = "\t")

ggplot(jSplice25, aes(x = logBaseMean ) ) + geom_histogram() 

## Create vector of geneIDs and see overlap with other programs
jSpliceVec <- jSplice25$GeneID

## Check overlap with DEGs
intersect(DEGs$ID, jSplice25$GeneID)
# [1] "AT3G13570" (SCL30a)

sum(jSpliceVec %in% finalSpladderGeneExpCounts$GeneID)
# 1 in common with spladder
jSpliceVec[jSpliceVec %in% finalSpladderGeneExpCounts$GeneID]
# [1] "AT1G61240"

## rMATS
sum(jSpliceVec %in% finalRMatsGeneExpCounts$GeneID)
# 12 in common with rMATs
jSpliceVec[jSpliceVec %in% finalRMatsGeneExpCounts$GeneID]
# [1] "AT4G39030" "AT3G59430" "AT4G22350" "AT3G13570" "AT2G17730" "AT1G11960" "AT5G65685" "AT3G22990" "AT5G50780"
# [10] "AT3G27870" "AT2G20950" "AT1G13990"

## SUPPA
sum(jSpliceVec %in% finalSuppaGeneExpCounts$GeneID)
# 20
jSpliceVec[jSpliceVec %in% finalSuppaGeneExpCounts$GeneID]
# [1] "AT1G68690" "AT4G39030" "AT3G59430" "AT4G01550" "AT1G14590" "AT4G30570" "AT4G32480" "AT1G11960" "AT5G24910"
# [10] "AT1G61240" "AT1G76570" "AT3G22990" "AT4G39350" "AT5G61830" "AT1G63670" "AT4G29000" "AT2G41720" "AT2G20950"
# [19] "AT1G44446" "AT3G26730"

sharedJspliceAndRmatsandSuppa <- intersect(jSpliceVec[jSpliceVec %in% finalRMatsGeneExpCounts$GeneID], 
          jSpliceVec[jSpliceVec %in% finalSuppaGeneExpCounts$GeneID])

# [1] "AT4G39030" "AT3G59430" "AT1G11960" "AT3G22990" "AT2G20950"

##############################################################
# FILTERING SPLADDER RESULTS ACCORDING TO jSPLICE PARAMETERS #
##############################################################

### Filtering spladder table based on min log2baseMean of 2.0 (equivalent to lowest baseMean in jSplice)
stringentSpladder <- finalSpladderGeneExpCounts %>% filter(logBaseMean > 2)

## Let's filter retained introns only and apply similar jSplice count threshold (20)
## intron_conf: number of spliced alignments spanning the intron
## will need to filter so that I catch the cases where RI is in the WT but not in KO and vice
## versa

RI_stringentSpladder <- stringentSpladder %>% filter(Type == "intron_retention") %>% 
        select_if(colSums(!is.na(.)) > 0) %>% 
        filter(gene_name != "AT2G08986", gene_name != "AT2G07981",
                                                 gene_name != "AT5G02500")

RI_WT_events <- RI_stringentSpladder %>% filter(atRWT2S.intron_conf > 20, atRWT3S.intron_conf > 20)
RI_KO_events <- RI_stringentSpladder %>% filter(atRKO1S.intron_conf > 20, atRKO3S.intron_conf > 20)

## Filter exon skipping cases, exon_pre_exon_aft_conf: number of spliced alignments spanning from left
## flanking to right flanking exon. 

## PROBLEM! The exon_pr_exon_aft_conf variable does not exist in my data. FRUSTRATING INCONCISTENCY
## WITH THE SPLADDER WIKI AND THE OUTPUT DATA. However, after closer inspection in IGV, intron_skip_conf is 
## equivalent to junction spanning reads that support skipping of the exon in question. AT4G13460 is
## a good example.

## remove NA only columns as well
ES_stringentSpladder <- stringentSpladder %>% filter(Type == "exon_skip") %>% 
        select_if(colSums(!is.na(.)) > 0)

ES_WT_events <- ES_stringentSpladder %>% filter(atRWT2S.intron_skip_conf > 20, atRWT3S.intron_skip_conf > 20)
ES_KO_events <- ES_stringentSpladder %>% filter(atRKO1S.intron_skip_conf > 20, atRKO3S.intron_skip_conf > 20)

## remove artefactual events based on crazy read alignments in IGV, i.e. reads that span huge distances
## Some of these reads span > 500,000 bases. HUGE insert sizes. Why?
## Furthermore, HSC70-1 has almost 120k reads mapped to it. Will exclude.
ES_KO_events <- ES_KO_events %>% filter(gene_name != "AT2G08986", gene_name != "AT2G07981",
                                        gene_name != "AT5G02500")
ES_WT_events <- ES_WT_events %>% filter(gene_name != "AT2G08986", gene_name != "AT2G07981",
                                        gene_name != "AT5G02500")

outersect <- function(x, y) {
        sort(c(setdiff(x, y),
               setdiff(y, x)))
}

## Get genes in WT ES events not in KO ES events
outersect(ES_KO_events$gene_name, ES_WT_events$gene_name)
setdiff(ES_WT_events$gene_name, ES_KO_events$gene_name )

## Now to filter the A5 and A3 events
A5_stringentSpladder <- stringentSpladder %>% filter(Type == "alt_5prime", gene_name != "AT2G08986", 
                                                     gene_name != "AT2G07981",
                                                     gene_name != "AT5G02500") %>% 
        select_if(colSums(!is.na(.)) > 0)

A3_stringentSpladder <- stringentSpladder %>% filter(Type == "alt_3prime", gene_name != "AT2G08986", 
                                                     gene_name != "AT2G07981",
                                                     gene_name != "AT5G02500") %>% 
        select_if(colSums(!is.na(.)) > 0)

A5_WT_events <- A5_stringentSpladder %>% 
        filter(atRWT2S.intron1_conf > 20, atRWT3S.intron2_conf > 20)

A5_KO_events <- A5_stringentSpladder %>% 
        filter(atRKO1S.intron1_conf > 20, atRKO3S.intron2_conf > 20)

A3_WT_events <- A3_stringentSpladder %>% 
        filter(atRWT2S.intron1_conf > 20, atRWT3S.intron2_conf > 20)

A3_KO_events <- A3_stringentSpladder %>% 
        filter(atRKO1S.intron1_conf > 20, atRKO3S.intron2_conf > 20)

## Mult_ex skip

multES_stringentSpladder <- stringentSpladder %>% filter(Type == "mult_exon_skip", gene_name != "AT2G08986", 
                                                         gene_name != "AT2G07981",
                                                         gene_name != "AT5G02500") %>% 
        select_if(colSums(!is.na(.)) > 0)

multES_stringentSpladder %>% select(contains("_conf")) %>% select_if(colSums(.) > 20)

# No gene satisfying criteria
