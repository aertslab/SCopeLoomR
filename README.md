# SCopeLoomR v0.2.2
An R package (compatible with SCope) to create generic .loom files and extend them with other data e.g.: SCENIC regulons, Seurat clusters and markers, ...

## Installation

Installation should take less than one minute.

```
install.packages("devtools")
library(devtools)
install_github("aertslab/SCopeLoomR")
```

## Tutorial
You can find a tutorial on how to create .loom files [here](https://github.com/aertslab/SCopeLoomR/blob/master/vignettes/SCopeLoomR_tutorial.Rmd).

## Version History

October 31, 2018

* Version 0.2.2
    * Add sanity tests when adding embeddings and SCENIC regulons with threshold assignments.
    * Small bug fixes.

August 8, 2018

* Version 0.2.1
    * Fix encoding bug leading to annotation not shown in SCope. 
    * Add sanity test to check that all cells from the dgem are present in the Seurat clusters.

July 6, 2018

* Version 0.2.0
    * Add feature to store trajectory data. Currently, the following softwares can have their trajectory data easily stored in .loom files: `monocle`.

July 4, 2018

* Version 0.1.0
    * Add feature to store metrics
        
