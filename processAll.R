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
dcut <- 0.3

## Get entire list of gene IDs from ATRT2D in order to query thalemine for gene symbols
## cat AtRTD2_19April2016.gtf | gawk 'BEGIN {FS = ";"} {print $2}' | uniq | gawk '{print $2}' | 
## sed 's/\"//g' > thalemine.genelist.txt

## Create lists for applying various functions
# This doesn't work for some reason. Will have to use Sys.glob
#rmatsList <- list.files(pattern = paste0(dir_prefix,"atrt2d_MATS/MATS_output/*.MATS.ReadsOnTargetAndJunctionCounts.txt"), 
#                                        full.names = TRUE)

rmatsList <- Sys.glob(paste0(dir_prefix, "atrt2d_MATS/MATS_output/*.MATS.ReadsOnTargetAndJunctionCounts.txt"))
suppaList <- Sys.glob(paste0(dir_prefix, "suppa/atrtd2_quasi/*.dpsi"))
suppaPsiVecList <- Sys.glob(paste0(dir_prefix, "suppa/atrtd2_quasi/*.psivec"))
spladderList_mgraphs <- Sys.glob(paste0(dir_prefix, "spladder/*.confirmed.txt"))
spladderList_testR <- Sys.glob(paste0(dir_prefix, "spladder/*.tsv"))

## Run the various processing functions on each list

## RMATS PROCESSING ####################################################
myRmatsEvents <- lapply(rmatsList, processRMATS, fdr, dpsi, dcut)
myRmatsGenes <- lapply(rmatsList, processRMATS, fdr, dpsi, dcut, events = FALSE)
myRmatsCoordsList <- lapply(rmatsList, rMATsCoords, fdr, dpsi, dcut = 2)

## Create output table ranked by deltaPSI
## The following tsv file was obtained from Araport by uploading the gene IDs found in the AtRTD2 April 2016 GTF file
geneNames3 <- tbl_df(read.table("allthalemine.genes.txt.tsv", stringsAsFactors = FALSE, header = FALSE,
                                sep = "\t", na.strings = "", quote = "", comment.char = "")) %>%
        dplyr::rename(GeneID = V1, symbol = V2, description = V3, is_obsolete = V4)


paulaTableEventsR <- bind_rows(myRmatsEvents) %>% dplyr::arrange(desc(abs(IncLevelDifference)))
finalTableEventsR <- left_join(paulaTableEventsR, geneNames3, by = "GeneID")
write.table(finalTableEventsR, file = paste0("AllEvents_ranked_by_dPSI_rMATS_dcut_", dcut, ".txt"), quote = FALSE, sep = "\t", row.names = FALSE)


########################################################################
## Spladder processing

mySpladderEvents <- mapply(processSpladder, spladderList_mgraphs, spladderList_testR, dpsi = dpsi, dcut = dcut, SIMPLIFY = FALSE)
mySpladderGenes <- mapply(processSpladder, spladderList_mgraphs, spladderList_testR, dpsi = dpsi, dcut = dcut, events = FALSE, 
                          SIMPLIFY = FALSE)

mySpladderCoordList <- mapply(SpladderCoords, spladderList_mgraphs, spladderList_testR, dpsi = dpsi, dcut = 2, SIMPLIFY = FALSE)

## Create a table arranged by dPSI

paulaTableEventsSp <- bind_rows(mySpladderEvents) %>% dplyr::arrange(desc(abs(dPSI)))
finalTableEventsSp <- left_join(paulaTableEventsSp, geneNames3, by = "GeneID")
write.table(finalTableEventsSp, file = paste0("AllEvents_ranked_by_dPSI_spladder_dcut_", dcut, ".txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE)

########################################################################

## Suppa processing

mySuppaEvents <- mapply(processSuppa, suppaList, suppaPsiVecList, fdr, dpsi, dcut, SIMPLIFY = FALSE)
mySuppaGenes <- mapply(processSuppa, suppaList, suppaPsiVecList, fdr, dpsi, dcut, events = FALSE, SIMPLIFY = FALSE)

## Create a table arranged by dPSI

paulaTableEventsSu <- bind_rows(mySuppaEvents) %>% dplyr::arrange(desc(abs(dPSI)))
finalTableEventsSu <- left_join(paulaTableEventsSu, geneNames3, by = "GeneID")
write.table(finalTableEventsSu, file = paste0("AllEvents_ranked_by_dPSI_SUPPA_", dcut, ".txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE)

########################################################################

## Venn diagram of unique genes across programs
## Will generate at http://venndiagram.res.oicr.on.ca
## Need to create a dataframe with geneIDs as variables for each program

head(paulaTableEvents)

forVennRmats <- finalTableEventsR$GeneID
forVennSuppa <- finalTableEventsSu$GeneID        
forVennSpladder <- finalTableEventsSp$GeneID

venn.diagram(list(rMATS=forVennRmats, SUPPA=forVennSuppa, SplAdder=forVennSpladder), filename = "3way_genes_updated.png",
             imagetype = "png", main = "Overlapping genes across programs", sub = "All Event Types, heuritistic filter on reps",
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
write.table(onlyrMats_suppa_table, file = paste0("AllEvents_ranked_by_dPSI_only_rmats_suppa_dcut_", dcut, ".txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE)
write.table(onlySuppa_rmats_table, file = paste0("AllEvents_ranked_by_dPSI_only_suppa_rmats_dcut_", dcut, ".txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE)

## Sanity check, unique genes should be 223
length(unique(onlyrMats_suppa_table$GeneID))
length(unique(onlySuppa_rmats_table$GeneID))

## With accompanying table with only geneIDs, have to join by hand to make excel table of uneven rows

write.table(commonAcross3, file = paste0("commonAcross3_dcut_", dcut, ".txt"), quote = FALSE, row.names = FALSE)
write.table(rmats_suppa, file = paste0("rmats_suppa_dcut_", dcut, ".txt"), quote = FALSE, row.names = FALSE)
write.table(rmats_spladder, file = paste0("rmats_spladder_dcut_", dcut, ".txt"), quote = FALSE, row.names = FALSE)
write.table(suppa_spladder, file = paste0("suppa_spladder_dcut_", dcut, ".txt"), quote = FALSE, row.names = FALSE)

## Pull out the genes common across all 3 programs, with each program's output

common10.rMATS <- finalTableEventsR[finalTableEventsR$GeneID %in% commonAcross3,]
common10.spladder <- finalTableEventsSp[finalTableEventsSp$GeneID %in% commonAcross3,]
common10.suppa <- finalTableEventsSu[finalTableEventsSu$GeneID %in% commonAcross3,]

write.table(common10.rMATS, file = paste0("AllEvents_ranked_by_dPSI_common10.rmats_dcut_", dcut, ".txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE)

write.table(common10.spladder, file = paste0("AllEvents_ranked_by_dPSI_common10.spladder_dcut_", dcut, ".txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE)

write.table(common10.suppa, file = paste0("AllEvents_ranked_by_dPSI_common10.suppa_dcut_", dcut, ".txt"), quote = FALSE, sep = "\t", 
            row.names = FALSE)


## Create graphs of events
## Order for SplAdder is: A3, A5, SE, RI, MSE
Types <- c("A3", "A5", "SE", "RI", "MSE")
Counts <- sapply(mySpladderEvents, nrow)
uniqGeneCounts <- sapply(mySpladderGenes, nrow) 

eventsDf <- tbl_df(data.frame(Types, Counts)) %>% dplyr::rename(types = Types, counts = Counts)
eventsDf$variable <- "events"

genesDf <- tbl_df(data.frame(Types, uniqGeneCounts)) %>% dplyr::rename(types = Types, counts = uniqGeneCounts)
genesDf$variable <- "genes"

combinedDf <- bind_rows(eventsDf, genesDf)

## Create plot of AS event types
png(paste0("spladderASCounts_dcut_", dcut, ".png"))
ggplot(combinedDf, aes(x = types, y = counts, fill = variable)) + 
        geom_bar(stat = "identity", position = "dodge") + ggtitle(paste0("splAdder counts of AS event types (dPSI > 0.2)\n total counts of AS events: ",
sum(Counts), ", AS genes: ", sum(uniqGeneCounts))) +
        labs( x = "AS event types" ) + geom_text(aes(label=counts), position = position_dodge(width = 1), vjust = 0) +
        guides(fill=guide_legend(title=NULL)) + theme_grey()
dev.off()

## rMATS
## Order is: A3, A5, MXE, RI, SE
Types <- c("A3", "A5", "MXE", "RI", "SE")
Counts <- sapply(myRmatsEvents, nrow)
uniqGeneCounts <- sapply(myRmatsGenes, nrow)

eventsDf <- tbl_df(data.frame(Types, Counts)) %>% dplyr::rename(types = Types, counts = Counts)
eventsDf$variable <- "events"

genesDf <- tbl_df(data.frame(Types, uniqGeneCounts)) %>% dplyr::rename(types = Types, counts = uniqGeneCounts)
genesDf$variable <- "genes"

combinedDf <- bind_rows(eventsDf, genesDf)

## Create plot of AS event types
png(paste0("rMATSASCounts_dcut_", dcut, ".png"))
ggplot(combinedDf, aes(x = types, y = counts, fill = variable)) + 
        geom_bar(stat = "identity", position = "dodge") + ggtitle(paste0("rMATS counts of AS event types (dPSI > 0.2)\n total counts of AS events: ",
                                                                         sum(Counts), ", AS genes: ", sum(uniqGeneCounts))) +
        labs( x = "AS event types" ) + geom_text(aes(label=counts), position = position_dodge(width = 1), vjust = 0) +
        guides(fill=guide_legend(title=NULL)) + theme_grey()
dev.off()

## SUPPA
## Order is: A3, A5, MXE, RI, SE
Types <- c("A3", "A5", "MXE", "RI", "SE")
Counts <- sapply(mySuppaEvents, nrow)
uniqGeneCounts <- sapply(mySuppaGenes, nrow)

eventsDf <- tbl_df(data.frame(Types, Counts)) %>% dplyr::rename(types = Types, counts = Counts)
eventsDf$variable <- "events"

genesDf <- tbl_df(data.frame(Types, uniqGeneCounts)) %>% dplyr::rename(types = Types, counts = uniqGeneCounts)
genesDf$variable <- "genes"

combinedDf <- bind_rows(eventsDf, genesDf)

## Create plot of AS event types
png(paste0("SUPPAsASCounts_dcut_", dcut, ".png"))
ggplot(combinedDf, aes(x = types, y = counts, fill = variable)) + 
        geom_bar(stat = "identity", position = "dodge") + ggtitle(paste0("SUPPA counts of AS event types (dPSI > 0.2)\n total counts of AS events: ",
                                                                         sum(Counts), ", AS genes: ", sum(uniqGeneCounts))) +
        labs( x = "AS event types" ) + geom_text(aes(label=counts), position = position_dodge(width = 1), vjust = 0) +
        guides(fill=guide_legend(title=NULL)) + theme_grey()
dev.off()

## Look at distributions of differences between replicates.
## Looks like a good cutoff would be 0.3

png("histogram_diff_psi_reps.png")

par(mfrow = c(3,2))

h <- hist(finalTableEventsR$dWT, plot = FALSE)
h$density <- h$counts/sum(h$counts)*100
plot(h, freq = FALSE, ylab = "percent", xlab = "difference in PSI between reps", main = "Distribution of differences\n in WT PSI for rMATS events")
abline(v = 0.2)
#rug(h$counts)

h <- hist(finalTableEventsR$dKO, plot = FALSE)
h$density <- h$counts/sum(h$counts)*100
plot(h, freq = FALSE, ylab = "percent", xlab = "difference in PSI between reps", main = "Distribution of differences\n in KO PSI for rMATS events")
abline(v = 0.2)


h <- hist(finalTableEventsSp$dWT, plot = FALSE)
h$density <- h$counts/sum(h$counts)*100
plot(h, freq = FALSE, ylab = "percent", xlab = "difference in PSI between reps", main = "Distribution of differences\n in WT PSI for SplAdder events")
abline(v = 0.2)

h <- hist(finalTableEventsSp$dKO, plot = FALSE)
h$density <- h$counts/sum(h$counts)*100
plot(h, freq = FALSE, ylab = "percent", xlab = "difference in PSI between reps", main = "Distribution of differences\n in KO PSI for SplAdder events")
abline(v = 0.2)


h <- hist(finalTableEventsSu$dWT, plot = FALSE)
h$density <- h$counts/sum(h$counts)*100
plot(h, freq = FALSE, ylab = "percent", xlab = "difference in PSI between reps", main = "Distribution of differences\n in WT PSI for SUPPA events")
abline(v = 0.2)


h <- hist(finalTableEventsSu$dKO, plot = FALSE)
h$density <- h$counts/sum(h$counts)*100
plot(h, freq = FALSE, ylab = "percent", xlab = "difference in PSI between reps", main = "Distribution of differences\n in KO PSI for SUPPA events")
abline(v = 0.2)

dev.off()

par(mfrow=c(1,1))

## Plot the PSI values against each other
## rMATS
plot(finalTableEventsR$atRWT2S.psi, finalTableEventsR$atRWT3S.psi, col = rgb(0, 0, 0, 0.2), pch = 19)
abline(lm(finalTableEventsR$atRWT3S.psi~finalTableEventsR$atRWT2S.psi), col="red") # regression line (y~x) 
cor(finalTableEventsR$atRWT2S.psi, finalTableEventsR$atRWT3S.psi)

plot(finalTableEventsR$atRKO1S.psi, finalTableEventsR$atRKO3S.psi, col = rgb(0, 0, 0, 0.2), pch = 19)

## SplAdder
plot(finalTableEventsSp$atRWT2S.psi, finalTableEventsSp$atRWT3S.psi, col = rgb(0, 0, 0, 0.2), pch = 19)
plot(finalTableEventsSp$atRKO1S.psi, finalTableEventsSp$atRKO3S.psi, col = rgb(0, 0, 0, 0.2), pch = 19)
cor(finalTableEventsSp$atRWT2S.psi, finalTableEventsSp$atRWT3S.psi)
cor(finalTableEventsSp$atRKO1S.psi, finalTableEventsSp$atRKO3S.psi)

## SUPPA
plot(finalTableEventsSu$atRWT2S.psi, finalTableEventsSu$atRWT3S.psi, col = rgb(0, 0, 0, 0.2), pch = 19)
plot(finalTableEventsSu$atRKO1S.psi, finalTableEventsSu$atRKO3S.psi, col = rgb(0, 0, 0, 0.2), pch = 19)
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
spladderRI.df <- mySpladderCoordList[[4]] %>% dplyr::select(contig:gene_name, intron_start, intron_end) %>% dplyr::rename(
        chr = contig, start = intron_start, end = intron_end, GeneID = gene_name)

spladderA3_1.df <- mySpladderCoordList[[1]] %>% dplyr::select(contig:gene_name, exon_alt1_start, exon_alt1_end) %>% dplyr::rename(
        chr = contig, start = exon_alt1_start, end = exon_alt1_end, GeneID = gene_name)

spladderA3_2.df <- mySpladderCoordList[[1]] %>% dplyr::select(contig:gene_name, exon_alt2_start, exon_alt2_end) %>% dplyr::rename(
        chr = contig, start = exon_alt2_start, end = exon_alt2_end, GeneID = gene_name)

grSpladderRI <- makeGRangesFromDataFrame(spladderRI.df, keep.extra.columns = TRUE)

grSpladderA31 <- makeGRangesFromDataFrame(spladderA3_1.df, keep.extra.columns = TRUE)
grSpladderA32 <- makeGRangesFromDataFrame(spladderA3_2.df, keep.extra.columns = TRUE)

#########################
#       SUPPA           #
#########################

## Build the gr object from the Suppa Table
SuppaRI.df <- finalTableEventsSu[finalTableEventsSu$Type == "RI",] %>% dplyr::select(event, GeneID)

## Chromosome
SuppaRI.df$chr <- sapply(strsplit(SuppaRI.df$event, ":"), "[", 2)
## RI coordinates
tmp <- sapply(strsplit(SuppaRI.df$event, ":"), "[", 4)

SuppaRI.df$start <- sapply(strsplit(tmp, "-"), "[", 1)
SuppaRI.df$end <- sapply(strsplit(tmp, "-"), "[", 2)

SuppaRI.df$strand <- sapply(strsplit(SuppaRI.df$event, ":"), "[", 6)
head(SuppaRI.df)

grSuppaRI <- makeGRangesFromDataFrame(SuppaRI.df, keep.extra.columns = TRUE)

###########################
#       RMATS             #
###########################

rMATSRI.df <- myRmatsCoordsList[[4]] %>% dplyr::select(chr, GeneID, ID, strand, riExonStart_0base, riExonEnd) %>%
        dplyr::rename(start = riExonStart_0base, end = riExonEnd)

nrow(rMATSRI.df)
# 508

rMATSA3_1.df <- myRmatsCoordsList[[1]] %>% dplyr::select(chr, GeneID, ID, strand, longExonStart_0base, longExonEnd) %>%
        dplyr::rename(start = longExonStart_0base, end = longExonEnd)

rMATSA3_2.df <- myRmatsCoordsList[[1]] %>% dplyr::select(chr, GeneID, ID, strand, shortES, shortEE) %>%
        dplyr::rename(start = shortES, end = shortEE)


grRmatsRI <- makeGRangesFromDataFrame(rMATSRI.df, keep.extra.columns = TRUE)

grRmatsA3_1 <- makeGRangesFromDataFrame(rMATSA3_1.df, keep.extra.columns = TRUE)

## Find overlaps between Spladder and Suppa and rMATS Intron Retention ###

spladderAndSuppa_overlaps <- findOverlaps(grSpladderRI, grSuppaRI)

overlaps_spladder_perspective_vs_Suppa_RI <- subsetByOverlaps(grSpladderRI, grSuppaRI, minoverlap = 10)

overlaps_spladder_perspective_vs_rMATS_RI <- subsetByOverlaps(grSpladderRI, grRmatsRI, minoverlap = 10)
# no overlaps!

overlaps_suppa_perspective_vs_rMATS_RI <- subsetByOverlaps(grSuppaRI, grRmatsRI, minoverlap = 10)
# 126 ranges

## Find overlaps between Spladder and Suppa and rMATS A3 SS

overlaps_spladder_perspective_vs_rMATS_A3 <- subsetByOverlaps(grSpladderA31, grRmatsA3_1, minoverlap = 10)
# 4 ranges
subsetByOverlaps(grSpladderA32, grRmatsA3_1, minoverlap = 10)
# same 4 ranges



