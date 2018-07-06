# SCopeLoomR v0.2.0
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

July 6, 2018

* Version 0.2.0
    * Changes:
        * Add feature to store trajectory data. Currently, the following softwares can have their trajectory data easily stored in .loom files: `monocle`.

July 4, 2018

* Version 0.1.0
    * Changes:
        * Add feature to store metrics
        
