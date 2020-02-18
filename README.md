# SCopeLoomR v0.6.4
An R package (compatible with SCope) to create generic .loom files and extend them with other data e.g.: SCENIC regulons, Seurat clusters and markers, ... The package can also be used to read data from .loom files.

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
# install.packages("devtools")
devtools::install_github("aertslab/SCopeLoomR")
```

## Tutorial
You can find a tutorial on how to create .loom files and read data from them in the [package vignette](https://github.com/aertslab/SCopeLoomR/blob/master/vignettes/SCopeLoomR_tutorial.Rmd).

## Version History

February, 2020

* Version 0.6.4
    * Fix bug build_loom not working if dgem is a data.frame

* Version 0.6.3
    * Fix bug opening loom file v2
    * Fix bug using add_hierarchy results in broken loom file for SCope
    * Fix conditional statement in get_dspace

* Version 0.6.2
    * Fix bug datasets should not use scalar space
    * Fix bug global_meta_data_exists to work for loom v3 specs
    * Fix bug conditional statement to create attrs for loom v3 specs
    * Fix bug use json value in init_global_meta_data

* Version 0.6.1
    * Fix bug LOOM_SPEC_VERSION attribute does not exist for old Loom files generated with SCopeLoomR version < 0.6.0.

* Version 0.6.0
    * Add compatibility for Loom v3
    * Fix bug getting default embedding

January, 2020

* Version 0.5.1
    * Fix bug when add_seurat_clustering with annotation argument

June, 2019

* Version 0.5.0
    * Adding Seurat clustering results through `add_seurat_clustering` function works now also with Seurat v3 objects.

February, 2019

* Version 0.4.0
    * Embeddings and trajectories inferred from TI methods available within the dyno framework can be stored in .loom files and be displayed in SCope.

January, 2019

* Version 0.3.5
    * close_loom(): New function
    * open_loom(): Added argument "mode" and option to open as read-only (mode="r")
    
* Version 0.3.4
    * New functions to read SCENIC results (get_regulonsAuc, get_embeddings, get_regulons, get_cellAnnotation, get_clusterings_withName, get_regulonThresholds).
    * Updates to documentation.
    
January 11, 2019

* Version 0.3.3
    * Fix bug when adding Seurat clustering.

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
        
