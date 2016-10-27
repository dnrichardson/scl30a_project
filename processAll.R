## October 25, 2016
## R script for post-processing output files of rMATS, SplAdder and SUPPA.


## Load in the helper functions
source("processing_functions.R")

library(dplyr)
library(ggplot2)
library(VennDiagram)
library(plyr)
library(tidyr)

dir_prefix <- "~/Dropbox/Current/newscl30a/"

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
myRmatsEvents <- lapply(rmatsList, processRMATS, 0.05, 0.2)
myRmatsGenes <- lapply(rmatsList, processRMATS, 0.05, 0.2, events = FALSE)

## Order is: A3, A5, MXE, RI, SE
Types <- c("A3", "A5", "MXE", "RI", "SE")
Counts <- sapply(myRmatsEvents, nrow)
uniqGeneCounts <- sapply(myRmatsGenes, nrow)

eventsDf <- tbl_df(data.frame(Types, Counts)) %>% rename(types = Types, counts = Counts)
eventsDf$variable <- "events"

genesDf <- tbl_df(data.frame(Types, uniqGeneCounts)) %>% rename(types = Types, counts = uniqGeneCounts)
genesDf$variable <- "genes"

combinedDf <- bind_rows(eventsDf, genesDf)

## Create output table ranked by deltaPSI
## The following tsv file was obtained from Araport by uploading the gene IDs found in the AtRTD2 April 2016 GTF file
geneNames3 <- tbl_df(read.table("allthalemine.genes.txt.tsv", stringsAsFactors = FALSE, header = FALSE,
                                sep = "\t", na.strings = "", quote = "", comment.char = "")) %>%
        dplyr::rename(GeneID = V1, symbol = V2, description = V3, is_obsolete = V4)


paulaTableEventsR <- bind_rows(myRmatsEvents) %>% dplyr::arrange(desc(abs(IncLevelDifference)))
finalTableEventsR <- left_join(paulaTableEventsR, geneNames3, by = "GeneID")
write.table(finalTableEventsR, file = "AllEvents_ranked_by_dPSI_rMATS.txt", quote = FALSE, sep = "\t", row.names = FALSE)


########################################################################
## Spladder processing

mySpladderEvents <- mapply(processSpladder, spladderList_mgraphs, spladderList_testR, dpsi = 0.2, SIMPLIFY = FALSE)
mySpladderGenes <- mapply(processSpladder, spladderList_mgraphs, spladderList_testR, dpsi = 0.2, events = FALSE, 
                          SIMPLIFY = FALSE)

## Create a table arranged by dPSI
#geneNames3 <- tbl_df(read.table("allthalemine.genes.txt.tsv", stringsAsFactors = FALSE, header = FALSE,
#                                sep = "\t", na.strings = "", quote = "", comment.char = "")) %>%
#        dplyr::rename(gene = V1, symbol = V2, description = V3, is_obsolete = V4)


paulaTableEventsSp <- bind_rows(mySpladderEvents) %>% dplyr::arrange(desc(abs(dPSI)))
finalTableEventsSp <- left_join(paulaTableEventsSp, geneNames3, by = "GeneID")
write.table(finalTableEventsSp, file = "AllEvents_ranked_by_dPSI_spladder.txt", quote = FALSE, sep = "\t", 
            row.names = FALSE)

########################################################################

## Suppa processing

mySuppaEvents <- mapply(processSuppa, suppaList, suppaPsiVecList, 0.05, 0.2, SIMPLIFY = FALSE)
mySuppaGenes <- mapply(processSuppa, suppaList, suppaPsiVecList, 0.05, 0.2, events = FALSE, SIMPLIFY = FALSE)

## Create a table arranged by dPSI

#geneNames3 <- tbl_df(read.table("allthalemine.genes.txt.tsv", stringsAsFactors = FALSE, header = FALSE,
#                                sep = "\t", na.strings = "", quote = "", comment.char = "")) %>%
#        dplyr::rename(geneID = V1, symbol = V2, description = V3, is_obsolete = V4)

paulaTableEventsSu <- bind_rows(mySuppaEvents) %>% dplyr::arrange(desc(abs(dPSI)))
finalTableEventsSu <- left_join(paulaTableEventsSu, geneNames3, by = "GeneID")
write.table(finalTableEventsSu, file = "AllEvents_ranked_by_dPSI_SUPPA.txt", quote = FALSE, sep = "\t", 
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
write.table(onlyrMats_suppa_table, file = "AllEvents_ranked_by_dPSI_only_rmats_suppa.txt", quote = FALSE, sep = "\t", 
            row.names = FALSE)
write.table(onlySuppa_rmats_table, file = "AllEvents_ranked_by_dPSI_only_suppa_rmats.txt", quote = FALSE, sep = "\t", 
            row.names = FALSE)

## Sanity check, unique genes should be 223
length(unique(onlyrMats_suppa_table$GeneID))
length(unique(onlySuppa_rmats_table$GeneID))

## With accompanying table with only geneIDs, have to join by hand to make excel table of uneven rows

write.table(commonAcross3, file = "commonAcross3.txt", quote = FALSE, row.names = FALSE)
write.table(rmats_suppa, file = "rmats_suppa.txt", quote = FALSE, row.names = FALSE)
write.table(rmats_spladder, file = "rmats_spladder.txt", quote = FALSE, row.names = FALSE)
write.table(suppa_spladder, file = "suppa_spladder.txt", quote = FALSE, row.names = FALSE)

## Pull out the genes common across all 3 programs, with each program's output

common10.rMATS <- finalTableEventsR[finalTableEventsR$GeneID %in% commonAcross3,]
common10.spladder <- finalTableEventsSp[finalTableEventsSp$GeneID %in% commonAcross3,]
common10.suppa <- finalTableEventsSu[finalTableEventsSu$GeneID %in% commonAcross3,]

write.table(common10.rMATS, file = "AllEvents_ranked_by_dPSI_common10.rmats.txt", quote = FALSE, sep = "\t", 
            row.names = FALSE)

write.table(common10.spladder, file = "AllEvents_ranked_by_dPSI_common10.spladder.txt", quote = FALSE, sep = "\t", 
            row.names = FALSE)

write.table(common10.suppa, file = "AllEvents_ranked_by_dPSI_common10.suppa.txt", quote = FALSE, sep = "\t", 
            row.names = FALSE)

