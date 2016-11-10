## Library of functions to process outputs of various differential AS testing programs.
## Spladder, rMATS and SUPPA


## Usage example:
## Create list of merge_graph files
# myList <- sort(list.files(pattern = "*.confirmed.txt", full.names = FALSE))

## Create list of test_results files
# testResults <- sort(list.files(pattern = ".tsv", full.names = FALSE))
## myEvents <- mapply(processSpladder, myList, testResults, dpsi = 0.2, SIMPLIFY = FALSE)
## myGenes <- mapply(processSpladder, myList, testResults, dpsi = 0.2, events = FALSE, SIMPLIFY = FALSE)

## The main function to process spladder output
processSpladder <- function (mergeGraphFile, testResultsFile, fdr = 0.05, dpsi = 0.1, dcut = 2, outfile = FALSE, events = TRUE,
                             dirs = TRUE){
        
        # Processes output from SplAdder and can generate human readable tables.
        #
        # Args:
        #       mergeGraphFile: The mergegraph output file of SplAdder
        #       testResultsFile: The *.tsv output file of spladder_test.py
        #       fdr: The false discovery rate cutoff you want to filter by
        #       dpsi: The delta PSI cutoff you want to filter by
        #       outfile: if TRUE, write output to file
        #       events: if TRUE, create output with all events, else condense events into unique gene IDs
        #       dirs: if TRUE, take the path into consideration when parsing out the "ASType"
        #               only use this if your input files to this function are not in the same working directory
        #       dcut: Default is 2, which returns all events regardless of concordance of PSI values in each replicate
        #       based on distributions of differences of psi values in my data, a cutoff of 0.3 should include most events
        
        ## Read in the merge_graph file, rename some variables and calculate dPSI
        confirmedEvents <- tbl_df(read.table(mergeGraphFile, stringsAsFactors = FALSE, header = TRUE)) %>% 
                setNames(gsub("Aligned.sortedByCoord.out","",names(.))) %>% 
                mutate(dPSI = 0.5*(atRWT2S.psi + atRWT3S.psi) - 0.5*(atRKO1S.psi + atRKO3S.psi))
        
        ## Read in the spladder_test.py results files
        dAS <- tbl_df(read.table(testResultsFile, stringsAsFactors = FALSE, header = TRUE)) %>% 
                dplyr::rename(GeneID = gene) %>%
                dplyr::select(event_id, GeneID, p_val_adj) %>%
                filter(p_val_adj < fdr)
                
        
        ## Take note of the AS type for file-naming purposes
        
        ## If there are directory paths involved, have to parse differently
        if (dirs == TRUE){
                
                tmp <- unlist(strsplit(mergeGraphFile, split = "\\." ))[1]
                
                ASType <- tail(unlist(strsplit(tmp, split = "/", fixed = TRUE)), n = 1)
                ASType <- gsub("merge_graphs_", "", ASType)
                ASType <- gsub("_C3", "", ASType)
        } else
                {
                ASType <- unlist(strsplit(mergeGraphFile, split = "\\." ))[1]
                ASType <- gsub("merge_graphs_", "", ASType)
                ASType <- gsub("_C3", "", ASType)
        }
        
        ## Join the matching rows based on event_id, filter with a heuristic to discard rows where difference
        ## between the replicates is greater than the delta psi. This is to ensure that our replicates agree
        
        eventswithdPSI <- inner_join(confirmedEvents, dAS, by = "event_id") %>%
                dplyr::mutate(Type = ASType, dWT = abs(atRWT2S.psi - atRWT3S.psi), dKO = abs(atRKO1S.psi - atRKO3S.psi)) %>%
                dplyr::select(event_id, GeneID, Type, atRWT2S.psi, atRWT3S.psi, atRKO1S.psi, atRKO3S.psi, dPSI, dWT, dKO, p_val_adj) %>% 
                dplyr::filter(abs(dPSI) > dpsi, abs(atRWT2S.psi - atRWT3S.psi) < dcut & abs(atRKO1S.psi - atRKO3S.psi) < dcut )
        
        
        if (events == TRUE){
                return (eventswithdPSI)
        } else if (events == FALSE) {
                ## Get a genewise list of unique ids
                uniqGenes <- distinct(eventswithdPSI, GeneID, .keep_all = TRUE)
                return (uniqGenes)
        }
        
        ## Write both to files
        if (outfile == TRUE & events == TRUE){
                write.table(eventswithdPSI, file = paste(ASType,"fdr", fdr, "dpsi", dpsi, "events.txt", sep = "_"),
                            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE )
        }
        else if (outfile == TRUE & events == FALSE) {
                write.table(uniqGenes,file = paste(ASType,"fdr", fdr, "dpsi", dpsi, "events.unique.gene.txt", sep = "_"),
                            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
        }
}

## Usage example
## Read in the delta psi files output by suppa.py diffSplice

# myList <- list.files(pattern = "*.dpsi", full.names = FALSE)

## myEvents <- lapply(myList, dPSIprocessor, 0.05, 0.2)
## myGenes <- lapply(myList, dPSIprocessor, 0.05, 0.2, events = FALSE)

## How many unique genes participating in AS?
# totalNumberofGeneswithAS <- bind_rows(myGenes) %>% distinct(geneID) %>% nrow

## Second spladder function just to get the coordinates of events

SpladderCoords <- function (mergeGraphFile, testResultsFile, fdr = 0.05, dpsi = 0.1, dcut = 2, outfile = FALSE,
                             dirs = TRUE){
        
        # Processes output from SplAdder and can generate human readable tables.
        #
        # Args:
        #       mergeGraphFile: The mergegraph output file of SplAdder
        #       testResultsFile: The *.tsv output file of spladder_test.py
        #       fdr: The false discovery rate cutoff you want to filter by
        #       dpsi: The delta PSI cutoff you want to filter by
        #       outfile: if TRUE, write output to file
        #       events: if TRUE, create output with all events, else condense events into unique gene IDs
        #       dirs: if TRUE, take the path into consideration when parsing out the "ASType"
        #               only use this if your input files to this function are not in the same working directory
        #       dcut: Default is 2, which returns all events regardless of concordance of PSI values in each replicate
        #       based on distributions of differences of psi values in my data, a cutoff of 0.3 should include most events
        
        ## Read in the merge_graph file, rename some variables and calculate dPSI
        confirmedEvents <- tbl_df(read.table(mergeGraphFile, stringsAsFactors = FALSE, header = TRUE)) %>% 
                setNames(gsub("Aligned.sortedByCoord.out","",names(.))) %>% 
                mutate(dPSI = 0.5*(atRWT2S.psi + atRWT3S.psi) - 0.5*(atRKO1S.psi + atRKO3S.psi))
        
        ## Read in the spladder_test.py results files
        dAS <- tbl_df(read.table(testResultsFile, stringsAsFactors = FALSE, header = TRUE)) %>% 
                dplyr::rename(GeneID = gene) %>%
                dplyr::select(event_id, GeneID, p_val_adj) %>%
                filter(p_val_adj < fdr)
        
        
        ## Take note of the AS type for file-naming purposes
        
        ## If there are directory paths involved, have to parse differently
        if (dirs == TRUE){
                
                tmp <- unlist(strsplit(mergeGraphFile, split = "\\." ))[1]
                
                ASType <- tail(unlist(strsplit(tmp, split = "/", fixed = TRUE)), n = 1)
                ASType <- gsub("merge_graphs_", "", ASType)
                ASType <- gsub("_C3", "", ASType)
        } else
        {
                ASType <- unlist(strsplit(mergeGraphFile, split = "\\." ))[1]
                ASType <- gsub("merge_graphs_", "", ASType)
                ASType <- gsub("_C3", "", ASType)
        }
        
        ## Join the matching rows based on event_id, filter with a heuristic to discard rows where difference
        ## between the replicates is greater than the delta psi. This is to ensure that our replicates agree
        
        ## Get coords
        
        eventCoords <- inner_join(confirmedEvents, dAS, by = "event_id") %>%
        dplyr::mutate(Type = ASType, dWT = abs(atRWT2S.psi - atRWT3S.psi), dKO = abs(atRKO1S.psi - atRKO3S.psi)) %>%
        dplyr::select(everything()) %>% 
        dplyr::filter(abs(dPSI) > dpsi, abs(atRWT2S.psi - atRWT3S.psi) < dcut & abs(atRKO1S.psi - atRKO3S.psi) < dcut )
        
        ## Get a genewise list of unique ids
        #uniqGenes <- distinct(eventswithdPSI, GeneID, .keep_all = TRUE)
        
        ## Write both to files
        if (outfile == TRUE){
                write.table(eventCoords, file = paste(ASType,"fdr", fdr, "dpsi", dpsi, "events_coords.txt", sep = "_"),
                            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE )
        }
        else{
                return (eventCoords)
}

}

## The main function to process suppa output
processSuppa <- function (dpsifile, psivecfile, fdr, dpsi = "", dcut = 2, outfile = FALSE, events = TRUE, dirs = TRUE){
        
        # Processes output from SUPPA and can generate human readable tables.
        #
        # Args:
        #       dpsifile: The *.dpsi output file of SUPPA diffSplice
        #       psivecfile: The *.psivec output file of SUPPA diffSplice
        #       fdr: The false discovery rate cutoff you want to filter by
        #       dpsi: The delta PSI cutoff you want to filter by
        #       outfile: if TRUE, write output to file
        #       events: if TRUE, create output with all events, else condense events into unique gene IDs
        #       dirs: if TRUE, take the path into consideration when parsing out the "ASType"
        #               only use this if your input files to this function are not in the same working directory
        
        ## Read in the delta psi files
        dPSIs <- tbl_df(read.table(dpsifile, skip = 1, stringsAsFactors = FALSE))
        
        ## Read in the psivec files to get individual psi values for each replicate
        
        psivecs <- tbl_df(read.table(psivecfile, stringsAsFactors = FALSE, header = FALSE, skip = 1)) %>% 
                dplyr::rename(event = V1, atRWT2S.psi = V2, atRWT3S.psi = V3, atRKO1S.psi = V4, atRKO3S.psi = V5)
        
        if (dirs == TRUE){
                ASType <- unlist(strsplit(dpsifile, split = "\\_" ))[3] %>% gsub(".dpsi", "",.)
        } else{
        
                ASType <- unlist(strsplit(dpsifile, split = "\\_" ))[2] %>% gsub(".dpsi", "",.)
        }
        
        # uniqIDs <- dPSIs %>%
        #         dplyr::filter(V3 < fdr & abs(V2) > dpsi ) %>%
        #         dplyr::mutate(geneID = sapply(strsplit(V1, ";"), "[", 1)) %>%
        #         dplyr::distinct(geneID, .keep_all = TRUE) %>%
        #         dplyr::rename(FDR = V3, dPSI = V2, event = V1, GeneID = geneID)
        
        allEvents <- dPSIs %>%
                dplyr::filter(V3 < fdr & abs(V2) > dpsi ) %>%
                dplyr::mutate(GeneID = sapply(strsplit(V1, ";"), "[", 1)) %>%
                dplyr::rename(FDR = V3, dPSI = V2, event = V1)
        
        allEvents$Type <- ASType
        
        ## Join the matching rows based on event_id
        eventswithdPSI <- inner_join(allEvents, psivecs, by = "event") %>%
                dplyr::mutate(dWT = abs(atRWT2S.psi - atRWT3S.psi), dKO = abs(atRKO1S.psi - atRKO3S.psi)) %>%
                dplyr::select(event, Type, GeneID, atRWT2S.psi, atRWT3S.psi, atRKO1S.psi, atRKO3S.psi, dPSI, dWT, dKO, FDR) %>%
                dplyr::filter(abs(dPSI) > dpsi, abs(atRWT2S.psi - atRWT3S.psi) < dcut & abs(atRKO1S.psi - atRKO3S.psi) < dcut )
        
        uniqIDs <- eventswithdPSI %>% distinct(GeneID, .keep_all = TRUE)
        
        if (outfile == TRUE){
                write.table(eventswithdPSI, file = paste(ASType,"fdr", fdr, "dpsi", dpsi, "events.txt", sep = "_"),
                            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE )
                
                write.table(uniqIDs, file = paste(myfile,"Uniq_atrtd_", fdr, "_", dpsi, ".txt"),
                            quote = FALSE, sep = "\t", row.names = FALSE )
        }
        else if (events == TRUE){
                return (eventswithdPSI)
        } else
                return (uniqIDs)
        
}

## This is the primary function that will take rMATS junction read count files and produce 
## condensed outputs that are more easy to work with

## Usage example:
## myEvents <- lapply(myList, sigEventprocessor, 0.05, 0.2)
## myGenes <- lapply(myList, sigEventprocessor, 0.05, 0.2, events = FALSE)

processRMATS <- function (myfile, fdr = 0.05, dpsi = "", dcut = 2, outfile = FALSE, events = TRUE, dirs = TRUE){
        
        # Processes output from rMATS and can generate human readable tables.
        #
        # Args:
        #       myfile: A "ReadsOnTargetAndJunctionCounts" output file of rMATS
        #       fdr: The false discovery rate cutoff you want to filter by
        #       dpsi: The delta PSI cutoff you want to filter by
        #       outfile: if TRUE, write output to file
        #       events: if TRUE, create output with all events, else condense events into unique gene IDs
        #       dirs: if TRUE, take the path into consideration when parsing out the "ASType"
        #               only use this if your input files to this function are not in the same working directory
        
        # TO DO: break apart inclusion levels into separate variables
        
        
        dPSIs <- tbl_df(read.table(myfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                                   quote = "\""))
       
        ## If there are directory paths involved, have to parse differently
        if (dirs == TRUE){
                tmp <- unlist(strsplit(myfile, split = "\\." ))[1]
                ASType <- unlist(strsplit(tmp, split = "/", fixed = TRUE))[9]
        } else{
                ASType <- unlist(strsplit(myfile, split = "\\." ))[1]
        }
        
        
        # uniqIDs <- dPSIs %>%
        #         filter(FDR < fdr & abs(IncLevelDifference) > dpsi ) %>%
        #         dplyr::select(GeneID, FDR:IncLevelDifference) %>%
        #         tidyr::separate(IncLevel1, into = c("atRWT2S.psi", "atRWT3S.psi"), sep = ",") %>%
        #         tidyr::separate(IncLevel2, into = c("atRKO1S.psi", "atRKO3S.psi"), sep = ",") %>%
        #         mutate_each_(funs(as.numeric), c("atRWT2S.psi", "atRWT3S.psi", "atRKO1S.psi", "atRKO3S.psi" )) %>%
        #         dplyr::filter(abs(IncLevelDifference) > abs(atRWT2S.psi - atRWT3S.psi) & abs(IncLevelDifference) > abs(atRKO1S.psi - atRKO3S.psi)) %>%
        #         distinct(GeneID, .keep_all = TRUE)
        
        allEvents <- dPSIs %>%
                filter(FDR < fdr & abs(IncLevelDifference) > dpsi ) %>% 
                dplyr::select(GeneID, ID, FDR:IncLevelDifference)
        
        allEvents$Type <- ASType
        
        allEvents <- dplyr::select(allEvents, GeneID, ID, Type, FDR:IncLevelDifference) %>% 
                tidyr::separate(IncLevel1, into = c("atRWT2S.psi", "atRWT3S.psi"), sep = ",") %>%
                tidyr::separate(IncLevel2, into = c("atRKO1S.psi", "atRKO3S.psi"), sep = ",") %>%
                mutate_each_(funs(as.numeric), c("atRWT2S.psi", "atRWT3S.psi", "atRKO1S.psi", "atRKO3S.psi" )) %>%
                dplyr::mutate(dWT = abs(atRWT2S.psi - atRWT3S.psi), dKO = abs(atRKO1S.psi - atRKO3S.psi)) %>%
                dplyr::select(GeneID, ID, Type, atRWT2S.psi, atRWT3S.psi, atRKO1S.psi, atRKO3S.psi, IncLevelDifference, dWT, dKO, FDR) %>%
                dplyr::filter(abs(atRWT2S.psi - atRWT3S.psi) < dcut & abs(atRKO1S.psi - atRKO3S.psi) < dcut )
        
        ## Get only unique Gene identifiers
        uniqIDs <- allEvents %>% distinct(GeneID, .keep_all = TRUE)
        
        if (outfile == TRUE){
                write.table(uniqIDs, file = paste(ASType,"fdr", fdr, "dpsi", dpsi, "uniq.genes.txt", sep = "_"),
                            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE )
                
                write.table(allEvents, file = paste(ASType,"fdr", fdr, "dpsi", dpsi, "all.events.txt", sep = "_"),
                            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE )
        } 
        else if (events == TRUE){
                return (allEvents)
        } else
                return (uniqIDs)
        
}

rMATsCoords <- function (myfile, fdr = 0.05, dpsi = "", dcut = 2, outfile = FALSE, events = TRUE, dirs = TRUE){
        
        # Processes output from rMATS and can generate human readable tables.
        #
        # Args:
        #       myfile: A "ReadsOnTargetAndJunctionCounts" output file of rMATS
        #       fdr: The false discovery rate cutoff you want to filter by
        #       dpsi: The delta PSI cutoff you want to filter by
        #       outfile: if TRUE, write output to file
        #       events: if TRUE, create output with all events, else condense events into unique gene IDs
        #       dirs: if TRUE, take the path into consideration when parsing out the "ASType"
        #               only use this if your input files to this function are not in the same working directory
        
        # TO DO: break apart inclusion levels into separate variables
        
        
        dPSIs <- tbl_df(read.table(myfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                                   quote = "\""))
        
        ## If there are directory paths involved, have to parse differently
        if (dirs == TRUE){
                tmp <- unlist(strsplit(myfile, split = "\\." ))[1]
                ASType <- unlist(strsplit(tmp, split = "/", fixed = TRUE))[9]
        } else{
                ASType <- unlist(strsplit(myfile, split = "\\." ))[1]
        }
        
        
        # uniqIDs <- dPSIs %>%
        #         filter(FDR < fdr & abs(IncLevelDifference) > dpsi ) %>%
        #         dplyr::select(GeneID, FDR:IncLevelDifference) %>%
        #         tidyr::separate(IncLevel1, into = c("atRWT2S.psi", "atRWT3S.psi"), sep = ",") %>%
        #         tidyr::separate(IncLevel2, into = c("atRKO1S.psi", "atRKO3S.psi"), sep = ",") %>%
        #         mutate_each_(funs(as.numeric), c("atRWT2S.psi", "atRWT3S.psi", "atRKO1S.psi", "atRKO3S.psi" )) %>%
        #         dplyr::filter(abs(IncLevelDifference) > abs(atRWT2S.psi - atRWT3S.psi) & abs(IncLevelDifference) > abs(atRKO1S.psi - atRKO3S.psi)) %>%
        #         distinct(GeneID, .keep_all = TRUE)
        
        allEvents <- dPSIs %>%
                filter(FDR < fdr & abs(IncLevelDifference) > dpsi ) %>% 
                dplyr::select(everything())
        
        allEvents$Type <- ASType
        
        allEvents <- allEvents %>% 
                tidyr::separate(IncLevel1, into = c("atRWT2S.psi", "atRWT3S.psi"), sep = ",") %>%
                tidyr::separate(IncLevel2, into = c("atRKO1S.psi", "atRKO3S.psi"), sep = ",") %>%
                mutate_each_(funs(as.numeric), c("atRWT2S.psi", "atRWT3S.psi", "atRKO1S.psi", "atRKO3S.psi" )) %>%
                dplyr::mutate(dWT = abs(atRWT2S.psi - atRWT3S.psi), dKO = abs(atRKO1S.psi - atRKO3S.psi)) %>%
                #dplyr::select(everything()) %>%
                dplyr::filter(abs(atRWT2S.psi - atRWT3S.psi) < dcut & abs(atRKO1S.psi - atRKO3S.psi) < dcut )
        
        if (outfile == TRUE){
        
                write.table(allEvents, file = paste(ASType,"fdr", fdr, "dpsi", dpsi, "all.events_coords.txt", sep = "_"),
                            quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE )
        } 
        else
                return (allEvents)
}


## This is a function that will create a volcano plot for rMATS output. It requires the calibrate library
## but right now it's a very prototypical function that needs some work to be generally useful.
makeVolcano <- function (inputFile, dpsi = 0.1, fdr = 0.05, annotate=FALSE) {
        png(paste(inputFile, "dpsi", dpsi,  "Volcano.png", sep = "_"))
        ASType <- unlist(strsplit(inputFile, split = "\\." ))[1]
        plotTitle <- paste(ASType, "dpsi > ", dpsi, "FDR < ", fdr)
        res <- tbl_df(read.table(inputFile, stringsAsFactors = FALSE, header = TRUE, 
                                 sep = "\t", quote = "\"")) 
        with(res, plot(IncLevelDifference, -log10(FDR), pch=20, main=plotTitle, xlim=c(-1.2,1.2)))
        with(subset(res, FDR<fdr & abs(IncLevelDifference)>dpsi), points(IncLevelDifference, -log10(FDR), pch=20, col="green"))
        
        if (annotate == TRUE){
                with(subset(res, FDR<fdr & abs(IncLevelDifference)>dpsi), textxy(IncLevelDifference, -log10(FDR), labs=substring(Event,1,9), cex=.4))
        }
        dev.off()
        
}

## Function to create a file for input to rmats2sashimi. This function will read in a junction read counts 
## file and output a file that contains output limited by FDR and dpsi 
makeSashimiInputs <- function(myfile, fdr = 0.05, dpsi = 0.1){
        dPSIs <- tbl_df(read.table(myfile, header = TRUE, sep = "\t", stringsAsFactors = FALSE,
                                   quote = "\"")) %>%
                filter(FDR < fdr, abs(IncLevelDifference) > dpsi)
        write.table(dPSIs, file = paste(myfile,"sashimi", fdr, dpsi, ".txt", sep = "_"), 
                    quote = FALSE, row.names = FALSE, 
                    col.names = TRUE, sep = "\t")
}

