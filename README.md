# SCopeLoomR v0.3.2
An R package (compatible with SCope) to create generic .loom files and extend them with other data e.g.: SCENIC regulons, Seurat clusters and markers, ... The package can also be used to extract
data from .loom files.

## Requirements
- HDF5 >= 1.10.1

For **Linux** and **MacOS** machines, version 1.10.1 can be installed with this snippet:
```
curl -O https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-1.10/hdf5-1.10.1/src/hdf5-1.10.1.tar
cd hdf5-1.10.1
./configure
make -j4
make check
make install
```
To install HDF5 on a **Windows** machine, please use the prebuilt binaries available at https://www.hdfgroup.org/downloads/hdf5/.


For other HDF5 releases, please visit https://support.hdfgroup.org/ftp/HDF5/releases/.

## Installation

Installation should take less than one minute.

```
install.packages("devtools")
library(devtools)
install_github("aertslab/SCopeLoomR")
```

## Tutorial
You can find a tutorial on how to create .loom files and extract data from them [here](https://github.com/aertslab/SCopeLoomR/blob/master/vignettes/SCopeLoomR_tutorial.Rmd).

## Version History

January 10, 2019

* Version 0.3.2
    * Fix bugs related to MetaData global attribute (absent or not compressed).
    * Add feature to overwrite default clustering if already existing one in the .loom.

November 8, 2018

* Version 0.3.1
    * Fix bug when adding sparse matrices.

* Version 0.3.0
    * Add feature to extract the gene expression matrix of a given cluster.
    * Update [tutorial](https://github.com/aertslab/SCopeLoomR/blob/master/vignettes/SCopeLoomR_tutorial.Rmd).

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
        
