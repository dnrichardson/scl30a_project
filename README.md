# Files for my ongoing research into SCL30a.

The file, "processing_functions.R" contains several functions that I have created and used to process output files arising from three alternative splicing programs:

* [SUPPA](https://bitbucket.org/regulatorygenomicsupf/suppa)
* [SplAdder](https://bitbucket.org/regulatorygenomicsupf/suppa)
* [rMATS](http://rnaseq-mats.sourceforge.net)

To use the functions, source the file into your current R session and be sure that you have dplyr loaded, i.e. 

```R
source("processing_functions.R")
library(dplyr)
```

There are currently 3 "usable" functions:

* processSpladder
* processSuppa
* processRMATS

Each are documented in the "processing_functions.R" file. An example of usage is the following:

```R
source("processing_functions.R")
library(dplyr)

## Create a list of input files for the processRMATS function
rmats_inputList <- Sys.glob("*.MATS.ReadsOnTargetAndJunctionCounts.txt")

## apply the function processRMATS to each file in the input list, setting FDR < 0.05 and Inclusion Level Difference > 0.2
myRmatsEvents <- lapply(rmats_inputList, processRMATS, 0.05, 0.2)

## apply the function processRMATS to each file in the input list, setting FDR < 0.05 and Inclusion Level Difference > 0.2
## but condense output to yield only unique genes
myRmatsGenes <- lapply(rmats_inputList, processRMATS, 0.05, 0.2, events = FALSE)


```
