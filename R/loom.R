# loom.R
# 

##############################
# Constants                  #
##############################

# Column attributes
CA_CELLID<-"CellID"
CA_DFLT_CLUSTERS_NAME<-"ClusterName"
CA_DFLT_CLUSTERS_ID<-"ClusterID"
CA_EMBEDDING_NAME<-"Embedding"
CA_EMBEDDING_DFLT_CNAMES<-c("_X","_Y")
CA_EXTRA_EMBEDDINGS_X_NAME<-"Embeddings_X"
CA_EXTRA_EMBEDDINGS_Y_NAME<-"Embeddings_Y"
CA_EXTRA_EMBEDDINGS_NAMES<-c(CA_EXTRA_EMBEDDINGS_X_NAME, CA_EXTRA_EMBEDDINGS_Y_NAME)
CA_CLUSTERINGS_NAME<-"Clusterings"

# Row attributes
RA_CLUSTERING_MARKERS_NAME<-"ClusterMarkers"
RA_GENE_NAME<-"Gene"

# Global attributes
GA_METADATA_NAME<-"MetaData"
GA_METADATA_ANNOTATIONS_NAME<-"annotations"
GA_METADATA_METRICS_NAME<-"metrics"
GA_METADATA_CLUSTERINGS_NAME<-"clusterings"
GA_METADATA_EMBEDDINGS_NAME<-"embeddings"
GA_METADATA_CLUSTERINGS_CLUSTER_MARKER_METRICS_NAME<-"clusterMarkerMetrics"
GA_TITLE_NAME<-"title"
GA_TITLE_GENOME<-"Genome"
GA_CREATION_DATE_NAME<-"CreationDate"
GA_R_VERSION_NAME<-"RVersion"

# Miscellaneous
BASE64_PATTERN<-"^([A-Za-z0-9+/]{4})*([A-Za-z0-9+/]{3}=|[A-Za-z0-9+/]{2}==)?$"

#*******************#
#*******************#
#  DATA INSERTION   #
#*******************#
#*******************#

###################################
# Global attribute data functions #
###################################

#'@title create_hierarchy
#'@description      Create hierarchy as named list of tree levels (SCopeTreeL1, SCopeTreeL2, SCopeTreeL3)
#'@param level.1.name  Name of top level of the hierarchy.
#'@param level.2.name  Name of the second level of the hierarchy.
#'@param level.3.name  Name of the third level of the hierarchy.
#'@export
create_hierarchy<-function(level.1.name = NULL
                         , level.2.name = NULL
                         , level.3.name = NULL) {
  if(is.null(level.1.name))
    stop("You need to define at least the first level when creating a hierarchy.")
  if(is.null(level.2.name) && !is.null(level.3.name))
    stop("You cannot define the third level while the second level is not defined.")
  if(is.null(level.1.name))
    level.1.name<-""
  if(is.null(level.2.name))
    level.2.name<-""
  if(is.null(level.3.name))
    level.3.name<-""
  return (list("SCopeTreeL1"=level.1.name, "SCopeTreeL2"=level.2.name, "SCopeTreeL3"=level.3.name))
}

#'@title add_hierarchy
#'@description      Add a 3-level hierarchy as global attributes to the given loom.
#'@param loom       The loom file handler.
#'@param hierarchy  A named list of the 3-levels names. Use create_hierarchy() to build it.
#'@export
add_hierarchy<-function(loom, hierarchy, overwrite = FALSE) {
  if(loom$mode=="r") stop("File open as read-only.")
  he<-hierarchy_exists(loom = loom)
  if(he & !overwrite) {
    stop("Hierarchy already exists for the given loom. You can overwrite the hierarchy in the given loom by the given hierarchy by setting overwrite option to TRUE.")
  }
  for(idx in seq_along(hierarchy)) {
    if(!he | !overwrite) {
      add_global_attr(loom = loom, key = names(hierarchy)[idx], value = hierarchy[[idx]])
    } else {
      update_global_attr(loom = loom, key = names(hierarchy)[idx], value = hierarchy[[idx]])
    }
  }
}

#'@title hierarchy_exists
#'@description  Check if hierarchy exists for the given loom.
#'@param loom   The loom file handler.
#'@export
hierarchy_exists<-function(loom) {
  return (loom$attr_exists(attr_name = "SCopeTreeL1") | loom$attr_exists(attr_name = "SCopeTreeL2") | loom$attr_exists(attr_name = "SCopeTreeL3"))
}

##############################
# Global Meta data functions #
##############################

#'@title add_global_md_metric
#'@description  Add the metric with the given name to the global MetaData attribute.
#'@param loom   The loom file handler.
#'@param name   The name of the annotation.
#'@param values The values of the annotation to be added.
add_global_md_metric<-function(loom
                             , name) {
  if(loom$mode=="r") stop("File open as read-only.")
  gmd<-get_global_meta_data(loom = loom)
  a<-gmd[[GA_METADATA_METRICS_NAME]]
  a[[length(a)+1]]<-list(name = name)
  gmd[[GA_METADATA_METRICS_NAME]]<-NULL
  gmd[[GA_METADATA_METRICS_NAME]]<-a
  update_global_meta_data(loom = loom, meta.data.json = rjson::toJSON(x = gmd))
}

#'@title add_global_md_annotation
#'@description  Add the annotation with the given name to the global MetaData attribute.
#'@param loom   The loom file handler.
#'@param name   The name of the annotation.
#'@param values The values of the annotation to be added.
add_global_md_annotation<-function(loom
                                 , name
                                 , values) {
  if(loom$mode=="r") stop("File open as read-only.")
  gmd<-get_global_meta_data(loom = loom)
  a<-gmd[[GA_METADATA_ANNOTATIONS_NAME]]
  a[[length(a)+1]]<-list(name = name, values = as.list(as.character(unique(values))))
  gmd[[GA_METADATA_ANNOTATIONS_NAME]]<-NULL
  gmd[[GA_METADATA_ANNOTATIONS_NAME]]<-a
  update_global_meta_data(loom = loom, meta.data.json = rjson::toJSON(x = gmd))
}

#'@title init_global_md_embeddings
#'@description  Remove all embeddings from the given loom.
#'@param loom   The loom file handler.
#'@export
init_global_md_embeddings<-function(loom) {
  gmd<-get_global_meta_data(loom = loom)
  gmd[[GA_METADATA_EMBEDDINGS_NAME]]<-NULL
  gmd[[GA_METADATA_EMBEDDINGS_NAME]]<-list()
  update_global_meta_data(loom = loom, meta.data.json = rjson::toJSON(x = gmd))
}

#'@title add_global_md_embedding
#'@description  Add the embedding with the given name to the global MetaData attribute.
#'@param loom       The loom file handler.
#'@param name       The name of the embedding to add.
#'@param is.default Is the given embedding the default embedding to use in the .loom file.
#'@param trajectory  A data.frame storing the X (first column) and Y (second column) coordinates of the trajectory.
add_global_md_embedding<-function(loom
                                , id
                                , name
                                , is.default = F
                                , trajectory = NULL) {
  if(loom$mode=="r") stop("File open as read-only.")
  gmd<-get_global_meta_data(loom = loom)
  e<-gmd[[GA_METADATA_EMBEDDINGS_NAME]]
  idx<-length(e)+1
  e[[idx]]<-list(id = as.character(id)
               , name = name
  )
  if(!is.null(trajectory)) {
    e[[idx]][["trajectory"]]<-trajectory
  }
  gmd[[GA_METADATA_EMBEDDINGS_NAME]]<-NULL
  gmd[[GA_METADATA_EMBEDDINGS_NAME]]<-e
  update_global_meta_data(loom = loom, meta.data.json = rjson::toJSON(x = gmd))
}

#'@title add_global_md_clustering_kv
#'@description  Add given value with the given key to the clustering with the given clustering.id of the global MetaData attribute
#'@param loom           The loom file handler.
#'@param clustering.id  The clustering ID of the clustering which the key/value is added to.
#'@param key            The name of the attribute to add.
#'@param value          The value of the attribute to add.
add_global_md_clustering_kv<-function(loom
                                    , clustering.id
                                    , key
                                    , value) {
  if(loom$mode=="r") stop("File open as read-only.")
  gmd<-get_global_meta_data(loom = loom)
  c<-gmd[[GA_METADATA_CLUSTERINGS_NAME]]
  
  # Filter the clusterings based on the given clustering.id
  mask<-lapply(c, function(x) {
    return (x[["id"]])
  }) == clustering.id
  
  if(sum(mask) == 0) {
    stop("Cannot add key/value to clustering that do not exists.")
  }
  idx<-which(mask == T)
  
  tmp.clustering<-c[mask][[1]]
  tmp.clustering[[key]]<-value
  gmd[[GA_METADATA_CLUSTERINGS_NAME]][[idx]]<-tmp.clustering
  update_global_meta_data(loom = loom, meta.data.json = rjson::toJSON(x = gmd))
  flush(loom = loom)
}

#'@title add_global_md_clustering
#'@description  Add the clustering annotation to the global MetaData attribute.
#'@param loom     The loom file handler.
#'@param group    The name of the group of clusterings.
#'@param name     The name given to the given clustering.
#'@param clusters A list of the the cluster id for each cell present in the same order as in the columns of gene expression matrix.
add_global_md_clustering<-function(loom
                                 , id 
                                 , group
                                 , name
                                 , clusters
                                 , annotation = NULL) {
  if(loom$mode=="r") stop("File open as read-only.")
  if(sum(is.na(names(clusters))) > 0) {
    stop("The names of the given clusters contains NAs.")
  }
  gmd<-get_global_meta_data(loom = loom)
  c<-gmd[[GA_METADATA_CLUSTERINGS_NAME]]
  if(is.factor(clusters)) {
    unique.clusters<-sort(as.integer(levels(clusters)), decreasing = F)
  } else {
    unique.clusters<-sort(unique(clusters), decreasing = F)
  }
  
  clusters<-lapply(X = unique.clusters, FUN = function(cluster.id) {
    description<-paste("NDA - Cluster", cluster.id)
    if(!is.null(annotation)) {
      # Force to have the same order
      annotation<-annotation[names(clusters)]
      # If annotation for the current cluster not empty then add
      # d<-annotation[annotation[[annotation.cluster.id.cl]] == cluster.id, annotation.cluster.description.cl]
      # Convert from factor to character vector to be used with nchar
      d<-as.character(unique(annotation[clusters == cluster.id])) 
      
      if(length(d) > 1) {
        stop("Annotation is not unique: multiple annotation correspond to a cluster ID.")
      }
      if(nchar(d)>0) {
        description<-paste0(d, " (",cluster.id,")")
      }
    }
    return (list(id = cluster.id
                 , description = description))
  })
  clusterings<-get_clusterings(loom = loom)
  clustering<-list(id = id
                 , group = group
                 , name = name
                 , clusters = clusters)
  c[[length(c)+1]]<-clustering
  gmd[[GA_METADATA_CLUSTERINGS_NAME]]<-NULL
  gmd[[GA_METADATA_CLUSTERINGS_NAME]]<-c
  update_global_meta_data(loom = loom, meta.data.json = rjson::toJSON(x = gmd))
  flush(loom = loom)
}

#'@title add_global_md_regulon_thresholds
#'@description  Add the regulon thresholds annotation to the global MetaData attribute.
#'@param loom                           The loom file handler.
#'@param regulon.threshold.assignments  The automated regulon cell assignments object generated by SCENIC.
#'@param regulon.enrichment.table       The regulon enrichment table generated by SCENIC.
add_global_md_regulon_thresholds <-function(loom
                                            , regulon.threshold.assignments
                                            , regulon.enrichment.table = NULL) {
  if(loom$mode=="r") stop("File open as read-only.")
  gmd<-get_global_meta_data(loom = loom)
  rT<-gmd[["regulonThresholds"]]
  if(!is.null(regulon.enrichment.table)) {
    # Store the motif name
    if(!("NES"%in%colnames(regulon.enrichment.table))) {
      stop("Error: NES column is not provided in the regulon enrichment table!")
    }
    if(!("motif"%in%colnames(regulon.enrichment.table))) {
      stop("Error: motif column is not provided in the regulon enrichment table!")
    }
    regulon.enrichment.table <- regulon.enrichment.table[with(regulon.enrichment.table, order(NES, decreasing = T)), ]
    regulon.enrichment.table <- regulon.enrichment.table[which(regulon.enrichment.table$TFinDB!=""),]
  }
  for(regulon in names(regulon.threshold.assignments)) {
    regulon.name<-strsplit(regulon, " ")[[1]][1]
    # Store the threshold assignments
    rta<-regulon.threshold.assignments[[regulon]]
    AUC.thresholds<-rta$aucThr$thresholds
    allThresholds<-do.call(what = "c", args = lapply(X = row.names(AUC.thresholds), FUN = function(threshold.name) {
      l<-list()
      l[[threshold.name]]<-AUC.thresholds[threshold.name,"threshold"]
      return (l)
    }))
    AUC.selected.threshold<-rta$aucThr$selected
    # Add regulon thresholds
    regulon.tresholds<-list(regulon = gsub(pattern = " ", replacement = "_", x = regulon)
                            , defaultThresholdValue = as.numeric(AUC.selected.threshold)
                            , defaultThresholdName = names(AUC.selected.threshold)
                            , allThresholds = allThresholds)
    # Add the motif data
    if(!is.null(regulon.enrichment.table)) {
      regulonTf <- gsub("_extended","",regulon.name)
      metSubset <- regulon.enrichment.table[regulon.enrichment.table$highlightedTFs == regulonTf,]
      motifByType <- split(metSubset$motif, metSubset$TFinDB)
      
      motif.name = NULL
      if(grepl("_extended",regulon.name))
      {
        motif.name = motifByType[["*"]][1] ## Extended
      }else{ 
        motif.name = motifByType[["**"]][1] ## High confidence annotation
      }
      
      motifData<-paste0(motif.name,".png")
      regulon.tresholds<-append(regulon.tresholds, list("motifData"=motifData))
    } else {
      warning("Argument 'regulon.enrichment.table' is not provided. Motif data will not be stored in the given .loom.")
    }
    rT[[length(rT)+1]]<-regulon.tresholds
  }
  gmd[["regulonThresholds"]]<-NULL
  gmd[["regulonThresholds"]]<-rT
  update_global_meta_data(loom = loom, meta.data.json = rjson::toJSON(x = gmd))
  flush(loom = loom)
}


#'@title check_global_meta_data
#'@description  Perform some sanity checks on the MetaData global attribute
#'@param loom   The loom file handler.
check_global_meta_data<-function(loom) {
  if(!global_meta_data_exists(loom = loom)) {
    stop("The MetaData global attribute does not exist in the given .loom. To be able to run this function please first initialize the MetaData global attribute by running init_global_meta_data(loom = loom).")
  }
}

#'@title global_meta_data_exists
#'@description  Check if global MetaData attribute exists for the given loom.
#'@param loom   The loom file handler.
global_meta_data_exists<-function(loom) {
  return (loom$attr_exists(attr_name = GA_METADATA_NAME))
}

#'@title get_global_meta_data
#'@description Get the global MetaData attribute as a R object.
#'@param loom The loom file handler.
#'@export
get_global_meta_data<-function(loom) {
  meta.data<-h5attr(x = loom, which = GA_METADATA_NAME)
  if(!is_base64_encoded(value = meta.data)) {
    if(!is_json(value = meta.data)) {
      stop("The global MetaData attribute in the given loom is corrupted.")
    }
  } else {
    meta.data<-decompress_gzb64(gzb64c = meta.data)
  }
  return (rjson::fromJSON(json_str = meta.data))
}

#'@title update_global_meta_data
#'@description Update the global MetaData attribute with the given meta.data.json.
#'@param loom           The loom file handler.
#'@param meta.data.json The meta data stored as a json string.
update_global_meta_data<-function(loom
                                  , meta.data.json) {
  compressed.meta.data<-compress_gzb64(c = as.character(meta.data.json))
  update_global_attr(loom = loom, key = GA_METADATA_NAME, value = as.character(compressed.meta.data))
}

#'@title init_global_meta_data
#'@description  Initialize the global MetaData attribute.
#'@param loom   The loom file handler.
#'@export
init_global_meta_data<-function(loom) {
  meta.data<-list(annotations = list()
                  , metrics = list()
                  , embeddings = list()
                  , clusterings = list()
                  , regulonThresholds = list())
  meta.data.json<-rjson::toJSON(meta.data)
  if(!(GA_METADATA_NAME %in% list.attributes(object = loom))) {
    compressed.meta.data<-compress_gzb64(c = as.character(meta.data.json))
    add_global_attr(loom = loom, key = GA_METADATA_NAME, value = as.character(compressed.meta.data))
  } else {
    stop("The MetaData global attribute already exists.")
  }
}

#'@title clear_global_meta_data
#'@description  Empty the global MetaData attribute.
#'@param loom   The loom file handler.
clear_global_meta_data<-function(loom) {
  meta.data<-list(annotations = list()
                , metrics = list()
                , embeddings = list()
                , clusterings = list()
                , regulonThresholds = list())
  meta.data.json<-rjson::toJSON(meta.data)
  if(!(GA_METADATA_NAME %in% list.attributes(object = loom))) {
    compressed.meta.data<-compress_gzb64(c = as.character(meta.data.json))
    add_global_attr(loom = loom, key = GA_METADATA_NAME, value = as.character(compressed.meta.data))
  } else {
    update_global_meta_data(loom = loom, meta.data.json = as.character(meta.data.json))
  }
}

#######################
# Embedding functions #
#######################

#'@title add_embedding_dyno
#'@description Add the given embedding as a column attribute and meta data related to the given embeddding to the given .loom file handler.
#'@param loom The loom file handler.
#'@param dyno A dyno object generated by dyno::infer_trajectory.
#'@param name A name to give to the embedding/trajectory inferred by the trajectory method.
#'@export
add_embedding_from_dyno<-function(loom, dyno, name = NULL) {
  if(is.null(x = name)) {
    stop("Please specify a name for the given analysis run with dyno.")
  }
  add_embedding(loom = loom
                , embedding = dyno$dimred
                , name = name
                , trajectory = create_trajectory_from_dyno(dyno = dyno)
                , is.default = !col_attrs_exists_by_key(loom = loom, key = "Embedding"))
}

#'@title add_default_embedding
#'@description  Add the main (default) embedding to the given loom file handler.
#'@param loom       The loom file handler.
#'@param embedding  The embedding to add.
add_default_embedding<-function(loom, embedding) {
  if(loom$mode=="r") stop("File open as read-only.")
  embedding<-as.data.frame(embedding)
  colnames(embedding)<-CA_EMBEDDING_DFLT_CNAMES
  add_col_attr(loom = loom, key = CA_EMBEDDING_NAME, value = embedding)
}

#'@title add_embedding
#'@description Add the given embedding as a column attribute and meta data related to the given embeddding to the given .loom file handler.
#'@param loom       The loom file handler.
#'@param embedding  A M-by-2 data.frame of the embeddings with M cells. 
#'                  If the given trajectory is not NULL, it should correspond to the projection of the cells onto the given trajectory.
#'@param name       The name of the given embedding.
#'@param is.default Default embedding to use in the .loom file.
#'@param trajectory A named list with the node names, the edges, and the x and y coordinates of the nodes of trajectory. Use the function create_trajectory().
#'@export
add_embedding<-function(loom
                        , embedding
                        , name
                        , is.default = F
                        , trajectory = NULL) {
  if(loom$mode=="r") stop("File open as read-only.")
  
  # Check if any NA values
  if(any(is.na(x = embedding))) {
    stop("Please make sure that the embeddings do not contain any NA values.")
  }
  
  # Check MetaData global attribute
  check_global_meta_data(loom = loom)
  
  if(is.null(row.names(embedding))) {
    stop("Please, make sure row.names are correctly filled in with cell IDs.")
  }
  
  # Check if some cells are missing betweem cells in embedding and the cells in the expression matrix
  check.missing.cells<-check_missings(x = get_cell_ids(loom = loom), y = row.names(x = embedding))
  if(!check.missing.cells$result) {
    if(length(check.missing.cells$x.missing.in.y) > 0) {
      stop(paste0("Some cells (",length(check.missing.cells$x.missing.in.y),") from expression matrix are missing in the embeddings: ", paste0(check.missing.cells$x.missing.in.y, collapse = ","), "."))
    }
    if(length(check.missing.cells$y.missing.in.x) > 0) {
      stop(paste0("Some cells (",length(check.missing.cells$y.missing.in.x),") from embeddings are missing in the expression matrix: ", paste0(check.missing.cells$y.missing.in.x, collapse = ","), "."))
    }
  }
  
  # Check if cells in embedding are in the same order as in the expression matrix
  if(!is_order_conserved(x = get_cell_ids(loom = loom), y = row.names(x = embedding))) {
    stop("Please, make sure that row.names of the given embedding are in the same order as the column names of the expression matrix in the given loom.")
  }
  
  # Add the default embedding also to Embeddings_X and Embeddings_Y
  if(is.default) {
    add_default_embedding(loom = loom, embedding = embedding)
  } else {
    print(paste0("Adding embedding ", name, "..."))
  }
  ca<-loom[["col_attrs"]]
  
  # Add a main embedding
  if(sum(CA_EXTRA_EMBEDDINGS_NAMES%in%names(ca)) != 2) {
    for(i in seq_along(CA_EXTRA_EMBEDDINGS_NAMES)) {
      e<-as.data.frame(embedding[,i])
      id<-"-1"
      colnames(e)<-id
      add_col_attr(loom = loom, key = CA_EXTRA_EMBEDDINGS_NAMES[i], value = e)
    }
  } else {
    for(i in seq_along(CA_EXTRA_EMBEDDINGS_NAMES)) {
      ca.embeddings<-get_col_attr_by_key(loom = loom, key = CA_EXTRA_EMBEDDINGS_NAMES[i])
      e<-as.data.frame(embedding[,i])
      id<-as.character(ncol(ca.embeddings))
      colnames(e)<-id
      ca.embeddings<-cbind(ca.embeddings, e)
      # Update the current coordinates Embeddings
      update_col_attr(loom = loom, key = CA_EXTRA_EMBEDDINGS_NAMES[i], value = as.data.frame(ca.embeddings))
    }
  }
  flush(loom = loom)
  # Add the given embedding to the global attribute MetaData
  add_global_md_embedding(loom = loom
                          , id = id
                          , name = name
                          , is.default = is.default
                          , trajectory = trajectory)
  flush(loom = loom)
}

##########################
# Trajectories functions #
##########################

#'@title create_trajectory_from_dyno
#'@description Create a trajectory object from dyno TI output.
#'@param dyno A dyno object generated by dyno::infer_trajectory.
#'@export
create_trajectory_from_dyno<-function(dyno) {
  # get cell positions
  cell_positions <- dyno$dimred %>% as.data.frame() %>% rownames_to_column("cell_id")
  
  # assign cells to closest milestone
  cell_positions <- left_join(
    cell_positions,
    dyno$milestone_percentages %>% group_by(cell_id) %>% arrange(desc(percentage)) %>% filter(dplyr::row_number() == 1) %>% select(cell_id, milestone_id),
    "cell_id"
  )
  
  print("Projecting the waypoints... (this may take some time)")
  waypoint_projection <- project_waypoints(
    traj = dyno,
    cell_positions = cell_positions
  )
  edges<-as.data.frame(x = waypoint_projection$edges[, c("from","to")]) # "arrow"
  coordinates<-data.frame("x"=waypoint_projection$positions[, "comp_1"],"y"=waypoint_projection$positions[, "comp_2"])
  # Set the row.names
  row.names(x = coordinates)<-waypoint_projection$positions$waypoint_id
  return (create_trajectory(edges = edges, coordinates = coordinates))
}

#'@title create_trajectory_from_monocle
#'@description Create a trajectory object from CellDataSet monocle output.
#'@param cds CellDataSet object from monocle.
#'@export
create_trajectory_from_monocle<-function(cds) {
  edges<-as.data.frame(x = igraph::get.edgelist(graph = cds@minSpanningTree))
  coordinates<-data.frame("x"=t(x = reducedDimK(cds = cds))[,1],"y"=t(x = reducedDimK(cds = cds))[,2])
  return (create_trajectory(edges = edges, coordinates = coordinates))
}

#'@title create_trajectory
#'@description Create a trajectory object from the given edges and coordinates of the nodes.
#'@param cds CellDataSet object from monocle.
#'@export
create_trajectory<-function(edges, coordinates) {
  if(!is.data.frame(x = edges))
    stop("The given edges should be a data.frame.")
  if(!is.data.frame(x = coordinates))
    stop("The given coordinates should be a data.frame.")
  if(all(as.character(seq(1,nrow(coordinates))) == row.names(x = coordinates)) | is.null(row.names(coordinates)))
    stop("The row.names of the given coordinates should have the node names.")
  colnames(edges)[c(1,2)]<-c("source","target")
  row.names(edges)<-NULL
  edges.list<-apply(edges, 1, as.list)
  nodes<-row.names(coordinates)
  colnames(coordinates)<-c("x","y")
  row.names(coordinates)<-NULL
  coord.list<-apply(coordinates, 1, as.list)
  return (list("nodes"=nodes
               , "edges"=edges.list
               , "coordinates"=coord.list))
}

#########################
# Clusterings functions #
#########################

#'@title get_seurat_clustering_resolutions
#'@description Get list of all computed clusterings resolutions in the given seurat object.
#'@param seurat A Seurat object.
#'@param prefix Prefix of the clustering to extract the resolutions from (only for Seurat v3 objects).
#'@export
get_seurat_clustering_resolutions<-function(seurat, prefix = NULL) {
  if(seurat@version >= 3) {
    clustering.ids <- grep(pattern = "res\\.", x = colnames(x = seurat@meta.data))
    clustering.prefix <- paste0(unique(x = stringr::str_split_fixed(string = colnames(x = seurat@meta.data)[clustering.ids], pattern = "res.", n = 2)[,1]), "res.")

    # Check if there are multiple prefixes for clustering information
    # Seurat v3 uses integrated_snn_ prefix for storing clustering information from its integration method
    if(length(x = clustering.prefix) > 1 & is.null(x = prefix)) {
      stop(paste0("Detected multiple prefixes for clustering information: ", paste(clustering.prefix, collapse = ", "), " . Please select one of them using the prefix argument."))
    } else if(length(x = clustering.prefix) == 1) {
      prefix <- clustering.prefix
    } else if(length(x = clustering.prefix) == 0) {
      stop("It seems that clustering data has not been computed. ")
    }
    if(length(x = prefix) == 1 & !(prefix %in% clustering.prefix)) {
      stop("The selected prefix does not exist.")
    }
    return (as.numeric(x = stringr::str_split_fixed(string = colnames(x = seurat@meta.data)[grep(pattern = paste0("^", prefix), x = colnames(x = seurat@meta.data))], pattern = prefix, n = 2)[,2]))
  } else if (seurat@version >= 2 & seurat@version < 3){
    return (as.numeric(x = stringr::str_split_fixed(string = names(seurat@calc.params)[grep(pattern = "FindClusters", x = names(seurat@calc.params))], pattern = "FindClusters.res.", n = 2)[,2]))
  } else {
    stop(paste0("This version of Seurat ", seurat@version," is not supported."))
  }
}

#'@title add_seurat_clustering
#'@description Add all the Seurat clusterings in the given seurat object to the given .loom file handler.
#'@param loom                                 The loom file handler.
#'@param seurat                               The Seurat object
#'@param seurat.assay                         The assay to access the data from (only if you're using Seurat version 3). Option used only for sanity checks.
#'@param seurat.markers.file.path.list        The named list of file paths to the markers saved in RDS format. The names should be the resolution id of the corresponding clustering (e.g.: res.2). Default is NULL. 
#'@param default.clustering.resolution        The clustering resolution (i.e.: res.2, ...) of the clustering that should be set as the default which an annotation can be set for.
#'@param annotation                           A data.frame with annotation for the clusters of the default clustering. Default is NULL.
#'@param annotation.cluster.id.cn             The column name to use for the IDs of the clusters found by the given clustering group. Default is NULL.
#'@param annotation.cluster.description.cn    The column name to use for the description of the clusters found by the given clustering group. Default is NULL.
#'@export
add_seurat_clustering<-function(loom
                                , seurat
                                , seurat.assay = "RNA"
                                , seurat.clustering.prefix = NULL
                                , seurat.markers.file.path.list = NULL
                                , seurat.marker.metric.accessors = NULL
                                , seurat.marker.metric.names = NULL
                                , seurat.marker.metric.description = NULL
                                , default.clustering.resolution = NULL
                                , default.clustering.overwrite = NULL
                                , annotation = NULL
                                , annotation.cluster.id.cn = NULL
                                , annotation.cluster.description.cn = NULL) {
  if(loom$mode=="r") stop("File open as read-only.")
  
  # Check MetaData global attribute
  check_global_meta_data(loom = loom)
  
  if(!is.null(seurat.markers.file.path.list)) {
    if(is.null(names(seurat.markers.file.path.list))) {
      stop("Argument 'seurat.markers.file.path.list' is not a named list. The names should correspond to the clustering ID in Seurat object (e.g.: res.2).")
    }
  }
  
  # Access data
  if(seurat@version >= 3) {
    counts <- seurat@assays[[seurat.assay]]@counts
    data <- seurat@assays[[seurat.assay]]@data
  } else if (seurat@version >= 2 & seurat@version < 3){
    counts <- seurat@raw.data
    data <- seurat@data
  } else {
    stop(paste0("This version of Seurat ", seurat@version," is not supported."))
  }
  
  # Check that all cells present in loom digital expression matrix are present in the given seurat@data and seurat@scale.data objects
  cells<-get_cell_ids(loom = loom)
  n.cells <- length(cells)
  if(sum(cells %in% colnames(x = counts)) != n.cells | sum(cells %in% colnames(x = data)) != n.cells) {
    stop(paste0("Some cells are missing. Please check that all cells from the digital expression matrix (dgem) in the given loom are present in the data slots i.e.: seurat@raw.data, seurat@data or seurat@assays[[", seurat.assay, "]]@counts, seurat@assays[[", seurat.assay, "]]@data"))
  }
  
  clustering.resolutions <- get_seurat_clustering_resolutions(seurat = seurat
                                                              , prefix = seurat.clustering.prefix)
  if(length(x = clustering.resolutions) > 0) {
    for(res in clustering.resolutions) {
      
      # Get cluster IDs for the current clustering
      resolution.id<-paste0(seurat.clustering.prefix,res)
      print(paste0("Seurat, ", seurat.clustering.prefix, res))
      
      if(seurat@version >= 3) {
        Seurat::Idents(object = seurat)<-paste0(seurat.clustering.prefix, res)
        cluster.ids<-Idents(object = seurat)
      } else if (seurat@version >= 2 & seurat@version < 3){
        seurat<-Seurat::SetAllIdent(object = seurat, id=resolution.id)
        cluster.ids<-seurat@ident
      } else {
        stop(paste0("This version of Seurat ", seurat@version," is not supported."))
      }
      
      is.default.clustering<-F
      
      # Add the Seurat clusters
      print("Adding Seurat clusters...")
      a<-NULL
      ac.id.cn<-NULL
      ac.description.cn<-NULL
      if(!is.null(x = default.clustering.resolution)) {
        if(res == default.clustering.resolution | resolution.id == default.clustering.resolution) {
          print("Adding default Seurat clusters...")
          if(!is.null(annotation) & !is.null(annotation.cluster.id.cn) & !is.null(annotation.cluster.description.cn)) {
            a<-annotation
            ac.id.cn<-annotation.cluster.id.cn
            ac.description.cn<-annotation.cluster.description.cn
          }
          is.default.clustering<-T
        }
      } else {
        # If only one clustering resolution computed then set this one as default
        if(length(clustering.resolutions) == 1) {
          is.default.clustering<-T
        }
      }
      flush(loom = loom)
      cluster.annotation<-create_cluster_annotation(clusters = cluster.ids, cluster.meta.data.df = a, cluster.id.cn = ac.id.cn,  cluster.description.cn = ac.description.cn)
      clid<-add_annotated_clustering(loom = loom
                                   , group = "Seurat"
                                   , name = paste("Seurat,", paste0(seurat.clustering.prefix, res))
                                   , clusters = cluster.ids
                                   , annotation = cluster.annotation
                                   , is.default = is.default.clustering
                                   , overwrite.default = default.clustering.overwrite)
      print(paste0("Clustering ID: ", clid))
      flush(loom = loom)
      # Add the Seurat markers if not empty
      if(!is.null(seurat.markers.file.path.list)) {
        print("Adding Seurat markers...")
        if(resolution.id %in% names(seurat.markers.file.path.list)) {
          if(file.exists(seurat.markers.file.path.list[[resolution.id]])) {
            seurat.markers<-readRDS(file = seurat.markers.file.path.list[[resolution.id]])
            seurat.markers.by.cluster<-split(x = seurat.markers, f = seurat.markers$cluster)
            add_clustering_markers(loom = loom
                                 , clustering.id = clid
                                 , clustering.markers = seurat.markers.by.cluster
                                 , marker.metric.accessors = seurat.marker.metric.accessors
                                 , marker.metric.names = seurat.marker.metric.names
                                 , marker.metric.descriptions = seurat.marker.metric.description)
            flush(loom = loom)
          } else {
            warning(paste0("Missing Seurat markers file for clustering resolution ", res, "."))
          }
        } else {
          warning(paste0("Seurat markers for clustering resolution ", res, " have not been computed."))
        }
      } else {
        print("No Seurat markers added.")
      }
    }
  } else {
    warning("It seems that you are using an old version of Seurat.")
    print("Adding Seurat clusters...")
    cluster.annotation<-create_cluster_annotation(clusters = seurat@ident, cluster.meta.data.df = annotation, cluster.id.cn = annotation.cluster.id.cn,  cluster.description.cn = annotation.cluster.description.cn)
    clid<-add_annotated_clustering(loom = loom
                                   , group = "Seurat"
                                   , name = "Seurat, resolution undefined"
                                   , clusters = seurat@ident
                                   , annotation = cluster.annotation
                                   , is.default = T
                                   , overwrite.default = default.clustering.overwrite)
  }
}

#'@title append_clustering_update_ca
#'@description  Append the given clustering with the given clustering.id.
#'@param loom           The loom file handler.
#'@param clustering.id  The id of the embedding group.
#'@param clustering     The embedding to add.
append_clustering_update_ca<-function(loom
                                      , clustering.id
                                      , clustering) {
  if(loom$mode=="r") stop("File open as read-only.")
  
  ca.clusterings<-get_clusterings(loom = loom)
  colnames(clustering)<-as.character(clustering.id)
  # Append this clustering
  ca.clusterings<-cbind(ca.clusterings, clustering)
  update_col_attr(loom = loom, key = CA_CLUSTERINGS_NAME, value = as.data.frame(x = ca.clusterings))
}

#'@title create_cluster_annotation
#'@description Create a named list where names are the cell IDs and the values are the annotation stored in the given cluster.description.cn column of the given cluster.meta.data.df. 
#'             If cluster.meta.data.df is set to NULL, annotation of cells are set to "NDA - Cluster X" where X is the cluster ID.
#'@param cluster.meta.data.df    A data.frame object with at least 2 columns named cluster.id.cn and cluster.description.cn.
#'@param clusters                A named list of the cell id and the assigned cluster id generated by a clustering method (e.g.: Seurat).
#'@param cluster.id.cn           The column name to use for the IDs of the clusters.
#'@param cluster.description.cn  The column name to use for the description of the clusters.
#'@export
create_cluster_annotation<-function(clusters
                                  , cluster.meta.data.df = NULL
                                  , cluster.id.cn =  NULL
                                  , cluster.description.cn = NULL) {
  if(is.factor(clusters)) {
    unique.clusters<-sort(as.integer(levels(clusters)), decreasing = F)
  } else {
    unique.clusters<-sort(unique(clusters), decreasing = F)
  }
  annotation<-setNames(object = rep(NA, length(clusters)), nm = names(clusters))
  for(cluster in unique.clusters) {
    description<-paste0("NDA - Cluster ", cluster)
    if(!is.null(cluster.meta.data.df)) {
      if(!(cluster.description.cn %in% colnames(cluster.meta.data.df))) {
        stop(paste0("The given column ",cluster.description.cn, " does not exists in the annotation provided."))
      }
      description<-cluster.meta.data.df[cluster.meta.data.df[[cluster.id.cn]] == cluster, cluster.description.cn]
    }
    annotation[clusters == cluster]<-description
  }
  annotation<-factor(x = annotation)
  names(x = annotation)<-names(clusters)
  return (annotation)
}

#'@title add_clustering
#'@description Add the given clusters in the given group column attribute and meta data related to the given clustering to the given .loom file handler.
#'@param loom               The loom file handler.
#'@param group              The for the given clustering group to which the given clusters have to be added.
#'@param name               The name given to this clustering.
#'@param clusters           A named list of the cell id and assigned the cluster id.
#'@param is.default         Set this clustering be set as default one.
#'@param overwrite.default  Overwrite the default clustering in the given loom if set to TRUE. Otherwise, it's not updated.
#'@export
add_clustering<-function(loom
                       , group
                       , name
                       , clusters
                       , is.default = F
                       , overwrite.default = NULL) {
  if(loom$mode=="r") stop("File open as read-only.")
  
  annotation<-create_cluster_annotation(clusters = clusters)
  add_annotated_clustering(loom = loom
                           , group = group
                           , name = name
                           , clusters = clusters
                           , annotation = annotation
                           , is.default = is.default
                           , overwrite.default = overwrite.default)
}

#'@title add_annotated_clustering
#'@description Add the given clusters in the given group column attribute and meta data related to the given clustering to the given .loom file handler.
#'@param loom                 The loom file handler.
#'@param group                The for the given clustering group to which the given clusters have to be added.
#'@param name                 The name given to this clustering
#'@param clusters             A named list of the cell id and assigned the cluster id.
#'@param annotation           A named list of the cell id and the corresponding annotation.
#'@param is.default           Set this clustering be set as default one.
#'@param overwrite.default    Overwrite the default clustering in the given loom if set to TRUE. Otherwise, it's not updated.
#'@export
add_annotated_clustering<-function(loom
                                 , group
                                 , name
                                 , clusters
                                 , annotation
                                 , is.default = F
                                 , overwrite.default = NULL) {
  if(loom$mode=="r") stop("File open as read-only.")
  
  id<-0
  if(length(unique(clusters)) == length(unique(annotation))) {
    # Make sure the order are the same
    annotation<-annotation[names(clusters)]
  } else {
    # Does not seem that cluster IDs and cluster annotation correspond
    # Remap to have the same number of unique IDs as the number of unique annotation
    library(plyr)
    clusters<-factor(x = as.integer(mapvalues(annotation, from = unique(x = annotation), to = seq_along(along.with = unique(x = annotation))-1)))
    names(clusters)<-names(annotation)
  }
  cell.ids<-get_cell_ids(loom = loom)
  # Check if all the cells are present in the given clusters
  n.mismatches<-sum(!(names(clusters) %in% cell.ids))
  if(n.mismatches > 0) {
    stop(paste0("Mismatches detected between the cell IDs (",n.mismatches,") in the given clusters object and in CellID column attribute of the .loom. Please do not use special characters for cell IDs except '_'."))
  }
  # Order the clusters in the order defined CellID column attribute
  clusters<-clusters[cell.ids]
  # If the clustering is the default one
  # Add it as the generic column attributes ClusterID and ClusterName
  if(is.default) {
    if(col_attrs_exists_by_key(loom = loom, key = CA_DFLT_CLUSTERS_ID) | col_attrs_exists_by_key(loom = loom, key = CA_DFLT_CLUSTERS_NAME)) {
      if(is.null(x = overwrite.default)) {
        stop("A default clustering has already been set. It can be overwritten by setting the corresponding *overwrite.default to TRUE or skipped by setting it to FALSE.")
      } else {
        if(overwrite.default) {
          warning("A default clustering has already been set. The current default clustering will be overwritten.")
          # Check for each if it exists. If it exists update otherwise create new.
          if(col_attrs_exists_by_key(loom = loom, key = CA_DFLT_CLUSTERS_ID)) {
            update_col_attr(loom = loom, key = CA_DFLT_CLUSTERS_ID, value = as.integer(as.character(x = clusters)))
          } else {
            add_col_attr(loom = loom, key = CA_DFLT_CLUSTERS_ID, value = as.integer(as.character(x = clusters)))
          }
          # Check if annotation is not null
          if(!is.null(x = annotation)) {
            if(col_attrs_exists_by_key(loom = loom, key = CA_DFLT_CLUSTERS_NAME))
              update_col_attr(loom = loom, key = CA_DFLT_CLUSTERS_NAME, value = as.character(x = annotation))
            else
              add_col_attr(loom = loom, key = CA_DFLT_CLUSTERS_NAME, value = as.integer(as.character(x = annotation)))
          } else {
            # As we overwrite ClusterID and annotation is null, if ClusterName is already present, to be consistent we should delete it.
            if(col_attrs_exists_by_key(loom = loom, key = CA_DFLT_CLUSTERS_NAME)) {
              remove_col_attr(loom = loom, key = CA_DFLT_CLUSTERS_NAME)
            }
          }
        } else {
          warning("A default clustering has already been set. It won't be overwritten.")
        }
      }
    } else {
      add_col_attr(loom = loom, key = CA_DFLT_CLUSTERS_ID, value = as.integer(as.character(x = clusters)))
      add_col_attr(loom = loom, key = CA_DFLT_CLUSTERS_NAME, value = as.character(x = annotation))
    }
  }
  # Adding the clustering data
  if(col_attrs_exists_by_key(loom = loom, key = CA_CLUSTERINGS_NAME)) {
    print(paste(CA_CLUSTERINGS_NAME, "already exists..."))
    ca.clusterings<-get_clusterings(loom = loom)
    # Set the clustering id
    id<-ncol(ca.clusterings) # n clusterings (start at 0)
    clustering<-data.frame("x" = as.integer(as.character(x = clusters)))
    append_clustering_update_ca(loom = loom, clustering.id = id, clustering = clustering)
  } else {
    print(paste(CA_CLUSTERINGS_NAME, "created..."))
    clustering<-data.frame("x" = as.integer(as.character(x = clusters)), stringsAsFactors = F)
    colnames(clustering)<-as.character(id)
    add_col_attr(loom = loom, key = CA_CLUSTERINGS_NAME, value = as.data.frame(x = clustering))
  }
  flush(loom = loom)
  
  # Adding the clustering meta data
  add_global_md_clustering(loom = loom
                           , id = id
                           , group = group
                           , name = name
                           , clusters = clusters
                           , annotation = annotation)
  flush(loom = loom)
  return (id)
}

####################
# SCENIC functions #
####################

#'@title add_scenic_regulons
#'@description Add the regulons with their target genes generated by SCENIC as a row attribute to the given .loom file handler.
#'@param loom                           The loom file handler.
#'@param regulons                       A list of list of the regulons and their target genes generated by SCENIC.
#'@param regulon.threshold.assignments  The automated regulon cell assignments object generated by SCENIC.
#'@param regulon.enrichment.table       The regulon enrichment table generated by SCENIC.
#'@export
add_scenic_regulons<-function(loom
                              , regulons
                              , regulon.threshold.assignments = NULL
                              , regulon.enrichment.table = NULL) {
  
  if(loom$mode=="r") stop("File open as read-only.")
  
  if(!is.null(regulon.threshold.assignments)) {
    cell.ids<-get_cell_ids(loom = loom)
    cells.not.found<-lapply(seq_along(along.with = regulon.threshold.assignments), function(i) {
      tf<-names(regulon.threshold.assignments)[i]
      x<-regulon.threshold.assignments[[i]]
      # Get all cells from regulon cell assignments not found in dgem
      cells.not.found<-x$assignment[!(x$assignment%in%cell.ids)]
      if(length(cells.not.found) > 0) {
        return (list(name=tf, cells.not.found=cells.not.found))
      }
    })
    # Remove all empty entries
    cells.not.found<-Filter(Negate(is.null), cells.not.found)
    if(length(cells.not.found) > 0) {
      msg<-sapply(X = cells.not.found, FUN = function(x) {
        return (paste0("Cells not found in expression matrix for ", x$name ," regulon: ", paste0(x$cells.not.found, collapse = ", ")))
      })
      cat(paste(msg, collapse="\n\n"))
    }
  }
  
  # Get cells IDs from expression matrix
  genes<-get_genes(loom = loom)
  
  # Add regulons
  regulons.mask<-do.call(what = "cbind", args = lapply(seq_along(regulons), function(regulon.idx) {
    reg.name<-names(regulons)[regulon.idx]
    reg.genes<-regulons[[regulon.idx]]
    reg.mask<-data.frame("x"=genes%in%reg.genes, stringsAsFactors = F)
    colnames(reg.mask)<-reg.name
    return (reg.mask)
  }))
  row.names(regulons.mask)<-genes
  colnames(regulons.mask)<-gsub(pattern = " ", replacement = "_", x = colnames(regulons.mask))
  add_row_attr(loom = loom, key = "Regulons", value = as.data.frame(x = regulons.mask))
  flush(loom = loom)
  
  # Add regulon thresholds
  if(!is.null(regulon.threshold.assignments)) {
    add_global_md_regulon_thresholds(loom = loom
                                   , regulon.threshold.assignments = regulon.threshold.assignments
                                   , regulon.enrichment.table = regulon.enrichment.table)
  }
  flush(loom = loom)
}

#'@title add_scenic_regulons_auc_matrix
#'@description Add the regulons AUC matrix generated by SCENIC as a column attribute to the given .loom file handler.
#'@param loom         The loom file handler.
#'@param regulons.AUC A matrix of the regulons AUC values with M regulons as rows and N cells as columns.
#'@export
add_scenic_regulons_auc_matrix<-function(loom
                                         , regulons.AUC) {
  
  if(loom$mode=="r") stop("File open as read-only.")
  
  # Check if some cells are missing betweem cells in embedding and the cells in the expression matrix
  check.missing.cells<-check_missings(x = get_cell_ids(loom = loom), y = colnames(x = regulons.AUC))
  if(!check.missing.cells$result) {
    if(length(check.missing.cells$x.missing.in.y) > 0) {
      stop(paste0("Some cells (",length(check.missing.cells$x.missing.in.y),") from expression matrix are missing in the regulon AUC matrix: ", paste0(check.missing.cells$x.missing.in.y, collapse = ","), "."))
    }
    if(length(check.missing.cells$y.missing.in.x) > 0) {
      stop(paste0("Some cells (",length(check.missing.cells$y.missing.in.x),") from regulon AUC matrix are missing in the expression matrix: ", paste0(check.missing.cells$y.missing.in.x, collapse = ","), "."))
    }
  }
  
  # Check if cells in embedding are in the same order as in the expression matrix
  if(!is_order_conserved(x = get_cell_ids(loom = loom), y = colnames(x = regulons.AUC))) {
    stop("Please, make sure that row.names of the given embedding are in the same order as the column names of the expression matrix in the given loom.")
  }
  add_col_attr(loom = loom, key = "RegulonsAUC", value = as.data.frame(x = t(regulons.AUC)))
  flush(loom = loom)
}

###########################
# Row data functions      #
###########################

#'@title add_clustering_markers
#'@description Add the clustering markers as a row attribute to the given .loom file handler.
#'@param loom                       The loom file handler.
#'@param dgem                       A matrix of the gene expression with M genes as rows and N cells as columns.
#'@param clustering.id              The clustering id that the given clustering.markers are specific for.
#'@param clustering.markers         A list of markers of each cluster found for the given clustering.id.
#'@param marker.metric.accessors    A list of the column names of the metrics to extract from the cluster marker data.frame in the given clustering.markers.
#'@param marker.metric.names        A list of names to attribute to the given marker.metric.accessor.
#'@param marker.metric.descriptions A list of description to attribute to the given marker.metric.accessor.
#'@export
add_clustering_markers<-function(loom
                               , clustering.id
                               , clustering.markers
                               , marker.metric.accessors = NULL
                               , marker.metric.names = NULL
                               , marker.metric.descriptions = NULL) {
  
  if(loom$mode=="r") stop("File open as read-only.")
  
  # Check if a gene column is present
  tmp<-colnames(clustering.markers[[1]])
  if(!("gene" %in% tmp)) {
    stop("Cannot find a 'gene' column. Please use 'gene' as column name")
  }
  # Add the clustering markers as a mask
  genes<-get_genes(loom = loom, is.flybase.gn = F)
  clustering.markers.mask<-do.call(what = "cbind", args = lapply(seq_along(clustering.markers), function(cluster.idx) {
    cluster.name<-names(clustering.markers)[cluster.idx]
    cluster.markers<-clustering.markers[[cluster.idx]][["gene"]]
    cm.mask<-data.frame("x" = genes %in% cluster.markers, stringsAsFactors = F)
    colnames(cm.mask)<-cluster.name
    return (cm.mask)
  }))
  row.names(clustering.markers.mask)<-genes
  print(paste0("Adding markers for clustering ", clustering.id, "..."))
  add_row_attr(loom = loom, key = paste0(RA_CLUSTERING_MARKERS_NAME, "_",clustering.id), value = as.data.frame(x = clustering.markers.mask))
  flush(loom = loom)
  
  # Add the marker metrics
  # Check all marker metric vectors have the same length
  if(!is.null(marker.metric.accessors) & !is.null(marker.metric.names) & !is.null(marker.metric.descriptions)) {
    if(!all(lapply(list(marker.metric.accessors,marker.metric.names,marker.metric.descriptions),length) == length(marker.metric.accessors))) {
      stop("The given marker.metric.accessors, marker.metric.names and marker.metric.descriptions should have the same length.")
    }
    print(paste0("Adding metrics for clustering ", clustering.id, "..."))
    metrics.av<-list()
    for(metric.idx in seq_along(along.with = marker.metric.names)) {
      metric.accessor<-marker.metric.accessors[metric.idx]
      metric.name<-marker.metric.names[metric.idx]
      metric.description<-marker.metric.descriptions[metric.idx]
      clustering.marker.metric<-do.call(what = "cbind", args = lapply(seq_along(clustering.markers), function(cluster.idx) {
        cluster.name<-names(clustering.markers)[cluster.idx]
        # Get the current metric.name in the cluster marker table of the current cluster.name
        cluster.markers<-clustering.markers[[cluster.idx]][, c("gene", metric.accessor)]
        genes.df<-data.frame("gene" = genes, stringsAsFactors = F)
        metric.df<-merge(x = genes.df, y = cluster.markers, by = "gene", all = T)
        metric.df[is.na(metric.df)] <- 0
        row.names(metric.df)<-metric.df$gene
        # Order by order of genes stored in .loom
        metric.df<-metric.df[match(genes, metric.df$gene),]
        metric.df[, "gene"]<-NULL
        colnames(metric.df)<-cluster.name
        return (metric.df)
      }))
      add_row_attr(loom = loom, key = paste0(RA_CLUSTERING_MARKERS_NAME, "_",clustering.id,"_",metric.accessor), value = as.data.frame(x = clustering.marker.metric))
      flush(loom = loom)
      metrics.av[[length(metrics.av)+1]]<-list("accessor"=metric.accessor, "name"=metric.name, "description"=metric.description)
    }
    add_global_md_clustering_kv(loom = loom, clustering.id = clustering.id, key = GA_METADATA_CLUSTERINGS_CLUSTER_MARKER_METRICS_NAME, value = metrics.av)
  }
}

#'@title add_fbgn
#'@description Add the Flybase gene as a row attribute to the given .loom file handler.
#'@param loom                       The loom file handler.
#'@param dgem                       A matrix of the gene expression with M genes as rows and N cells as columns.
#'@param fbgn.gn.mapping.file.path  A N-by-2 data.frame containing the mapping between the Flybase gene and the gene symbol.
#'@export
add_fbgn<-function(loom
                   , dgem
                   , fbgn.gn.mapping.file.path) {
  if(loom$mode=="r") stop("File open as read-only.")
  fbgn.gn.mapping<-utils::read.table(file = fbgn.gn.mapping.file.path, header = F, sep = "\t", quote = '', stringsAsFactors = F)
  colnames(fbgn.gn.mapping)<-c("FBgn",RA_GENE_NAME)
  tmp<-data.frame(row.names(dgem))
  colnames(tmp)<-RA_GENE_NAME
  genes<-merge(x = tmp, y = fbgn.gn.mapping, by = RA_GENE_NAME)
  add_row_attr(loom = loom, key = "FBgn", value = genes$FBgn)
}

#####################
# Generic functions #
#####################

#'@title lookup_all_global_attr
#'@description List all global attributes of the given .loom file handler.
#'@param loom The loom file handler.
#'@export
lookup_all_global_attr<-function(loom) {
  list.attributes(object = loom)
}


#'@title update_global_attr
#'@description Update the global attribute with the given key and the given value.
#'@param loom   The loom file handler.
#'@param key    The key of the global attribute to update.
#'@param value  The new value
update_global_attr<-function(loom
                             , key
                             , value) {
  if(loom$mode=="r") stop("File open as read-only.")
  remove_global_attr(loom = loom, key = key)
  add_global_attr(loom = loom, key = key, value = value)
}


#'@title remove_global_attr
#'@description Remove the global attribute with the given key.
#'@param loom The loom file handler.
#'@param key  The key of the global attribute to remove.
remove_global_attr<-function(loom
                             , key) {
  if(loom$mode=="r") stop("File open as read-only.")
  loom$attr_delete(attr_name = key)
  flush(loom = loom)
}


#'@title get_global_attr
#'@description Get the global attribute with the given key.
#'@param loom The loom file handler.
#'@param key  The key of the global attribute to get value from.
#'@export
get_global_attr<-function(loom
                          , key) {
  h5attr(x = loom, which = key)
}

#'@title add_global_attr
#'@description Add a new global attribute to the given .loom object accessible by the given key and containing the given value.
#'@param loom   The loom file handler.
#'@param key    The name of the new added attribute.
#'@param value  The value of the new added attribute.
add_global_attr<-function(loom
                          , key
                          , value
                          , dtype = NULL) {
  if(loom$mode=="r") stop("File open as read-only.")
  if(is.null(dtype)) {
    dtype<-guess_dtype(x = value)
  }
  # Encode character vector to UTF-8
  dtype<-hdf5_utf8_encode(value = value, dtype = dtype)
  loom$create_attr(attr_name = key, robj = value, dtype = dtype, space = get_dspace(x = "scalar"))
  flush(loom = loom)
}

#'@title remove_row_attr
#'@description  Remove the row attribute with the given key.
#'@param key    The name of the row attribute.
#'@export
remove_row_attr<-function(loom
                          , key) {
  if(loom$mode=="r") stop("File open as read-only.")
  loom$link_delete(name = paste0("row_attrs/", key))
  flush(loom = loom)
}

#'@title update_row_attr
#'@description  Update the row attribute with the given key and the given value.
#'@param key    The name of the row attribute to remove.
#'@param value  The new value to update.
update_row_attr<-function(loom
                          , key
                          , value) {
  if(loom$mode=="r") stop("File open as read-only.")
  loom$link_delete(name = paste0("row_attrs/",key))
  add_col_attr(loom = loom, key = key, value = value)
}

#'@title add_row_attr
#'@description Add a new row attribute to the given .loom object accessible by the given key and containing the given value.
#'@param loom   The loom file handler.
#'@param key    The name of the new added attribute.
#'@param value  The value of the new added attribute.
#'@export
add_row_attr<-function(loom
                       , key
                       , value
                       , dtype = NULL) {
  if(loom$mode=="r") stop("File open as read-only.")
  if(is.null(dtype)) {
    dtype<-guess_dtype(x = value)
  }
  # Encode character vector to UTF-8
  dtype<-hdf5_utf8_encode(value = value, dtype = dtype)
  loom$create_dataset(name = paste0("row_attrs/",key), robj = value, dtype = dtype)
  flush(loom = loom)
}

#'@title col_attrs_exists_by_key
#'@description Check if the column attribute with the given key exists.
#'@param loom The loom file handler.
#'@param key  The name of the new added attribute.
col_attrs_exists_by_key<-function(loom
                                , key) {
  ca<-loom[["col_attrs"]]
  return (key %in% names(ca))
}

#'@title get_col_attr_by_key
#'@description Get the column attribute with the given key.
#'@param loom The loom file handler.
#'@param key  The name of the column attribute.
#'@export
get_col_attr_by_key<-function(loom
                            , key) {
  ca<-loom[["col_attrs"]]
  return (ca[[key]][])
}

#'@title remove_col_attr
#'@description Remove the column attribute with the given key.
#'@param loom The loom file handler.
#'@param key  The name of column attribute to remove.
#'@export
remove_col_attr<-function(loom
                          , key) {
  if(loom$mode=="r") stop("File open as read-only.")
  loom$link_delete(name = paste0("col_attrs/", key))
  flush(loom = loom)
}

#'@title update_col_attr
#'@description Update the column attribute with the given key and with the given value.
#'@param loom   The loom file handler.
#'@param key    The name of column attribute to update.
#'@param value  The new value to update.
update_col_attr<-function(loom
                          , key
                          , value) {
  if(loom$mode=="r") stop("File open as read-only.")
  loom$link_delete(name = paste0("col_attrs/",key))
  add_col_attr(loom = loom, key = key, value = value)
}

#'@title add_col_attr
#'@description Add a new column attribute to the given .loom object accessible by the given key and containing the given value.
#'@param loom          The loom file handler.
#'@param key           The name of the new added attribute.
#'@param value         The value of the new added attribute.
#'@param as.annotation Define this attribute as a discrete attribute. This attribute will be visible in the compare tab.
#'@param as.metric     Define this attribute as a continues
#'@export
add_col_attr<-function(loom
                       , key
                       , value
                       , dtype = NULL
                       , as.annotation = F
                       , as.metric = F) {
  if(loom$mode=="r") stop("File open as read-only.")
  if(as.annotation & as.metric) {
    stop("Cannot add column attribute that is both of type metric and annotation.")
  }
  
  if(as.annotation & length(unique(value)) > 245) {
    stop("Cannot add column attribute as an annotation with more than 245 unique values.")
  }
  
  if(is.null(dtype)) {
    dtype<-guess_dtype(x = value)
  }
  # Encode character vector to UTF-8
  dtype<-hdf5_utf8_encode(value = value, dtype = dtype)
  loom$create_dataset(name = paste0("col_attrs/",key), robj = value, dtype = dtype)
  flush(loom = loom)
  if(as.annotation) {
    add_global_md_annotation(loom = loom, name = key, values = as.character(value))
    flush(loom = loom)
  }
  if(as.metric) {
    add_global_md_metric(loom = loom, name = key)
    flush(loom = loom)
  }
}

#'@title add_matrix
#'@description Add the given gene expression matrix dgem to the given .loom object.
#'@param loom The loom file handler.
#'@param dgem A matrix of the gene expression with M genes as rows and N cells as columns.
add_matrix<-function(loom
                     , dgem
                     , chunk.size
                     , display.progress) {
  if(loom$mode=="r") stop("File open as read-only.")
  
  # Check if any NA values
  if(any(is.na(x = dgem))) {
    stop("Please make sure that the expression matrix (dgem) does not contain any NA values.")
  }
  
  row.names(dgem)<-NULL
  colnames(dgem)<-NULL
  dtype<-get_dtype(x = dgem[1, 1])
  loom$create_dataset(
    name = 'matrix',
    dtype = dtype,
    dims = rev(x = dim(x = dgem))
  )
  chunk.points<-chunk_points(
    data.size = dim(x = dgem)[2],
    chunk.size = chunk.size
  )
  if (display.progress) {
    pb<-utils::txtProgressBar(char = '=', style = 3)
  }

  for (col in 1:ncol(x = chunk.points)) {
    row.start <- chunk.points[1, col]
    row.end <- chunk.points[2, col]
    loom[['matrix']][row.start:row.end, ] <- t(x = as.matrix(x = dgem[, row.start:row.end]))
    if(display.progress) {
      utils::setTxtProgressBar(pb = pb, value = col / ncol(x = chunk.points))
    }
  }
  flush(loom = loom)
}

#'@title flush
#'@description flush
#'@param loom The loom file handler.
#'@export
flush<-function(loom) {
  if(loom$mode=="r") stop("File open as read-only.")
  gc()
  loom$flush()
}

#'@title finalize
#'@description Save and close the given .loom file handler.
#'@param loom The loom file handler.
#'@export
finalize<-function(loom) {
  if(loom$mode=="r") {
    warning("File open as read-only. Only closing file (not flushed)...")
    flush(loom = loom)
  }
  loom$close_all()
}

#'@title close_loom
#'@description Close the given .loom file handler. Equivalent to finalize but without flushing.
#'@param loom The loom file handler.
#'@export
close_loom<-function(loom) {
  tryCatch(loom$close_all(), 
           error=function(e) {invisible(e)})
}

#'@title build_loom
#'
#'@description Write the data given as arguments to the given file.name .loom file.
#'@param file.name                  A string naming the .loom file to be generated.
#'@param title                      A short description of content of loom.
#'@param genome                     The genome used for the mapping.
#'@param dgem                       A matrix of the gene expression with M genes as rows and N cells as columns.
#'@param default.embedding          A M-by-2 data.frame of the embedding (X and Y coordinates) of the cells.
#'@param default.embedding.name     A description name for the given default.embedding
#'@param hierarchy                  A named list of the 3 hierarchy levels that can be used to group looms in SCope. Use create_hierarchy().
#'@param fbgn.gn.mapping.file.path  A N-by-2 data.frame containing the mapping between the Flybase gene and the gene symbol.
#'@param chunk.size                 The size of chunk of the gene expression matrix.
#'@param display.progress           Display progress when adding the gene expression matrix.
#'@export
build_loom<-function(file.name
                     , title = NULL
                     , genome = NULL
                     , dgem
                     , default.embedding = NULL
                     , default.embedding.name = NULL
                     , hierarchy = NULL
                     , fbgn.gn.mapping.file.path = NULL
                     , chunk.size = 1000
                     , display.progress = T) {
  is.dgem.sparse<-F
  loom<-H5File$new(filename = file.name, mode = "w")
  tryCatch({
    print("Adding global attributes...")
    # title
    if(!is.null(title)) {
      add_global_attr(loom = loom, key = GA_TITLE_NAME, value = as.character(title))
    }
    # Genome
    if(!is.null(genome)) {
      add_global_attr(loom = loom, key = GA_TITLE_GENOME, value = as.character(genome))
    }
    # Creation data
    add_global_attr(loom = loom, key = GA_CREATION_DATE_NAME, value = as.character(Sys.time()))
    # R version
    add_global_attr(loom = loom, key = GA_R_VERSION_NAME, value = as.character(R.version.string))
    
    cn<-colnames(dgem)
    rn<-row.names(dgem)
    # global MetaData attribute
    init_global_meta_data(loom = loom)
    
    # Add hierarchy levels
    if(!is.null(hierarchy)) {
      add_hierarchy(loom = loom, hierarchy = hierarchy)
    }
    # matrix
    # Check the type of the sparse matrix
    # convert to dgCMatrix if necessary to speedup populating the matrix slot
    if(class(dgem) == "dgTMatrix") {
      print("Converting to dgCMatrix...")
      dgem<-methods::as(object = dgem, Class = "dgCMatrix")
    }
    # Check if sparse matrix
    if(sum(class(x = dgem) %in% c("dgCMatrix","dgTMatrix")) > 0) {
      is.dgem.sparse<-T
    }
    print("Adding matrix...")
    add_matrix(loom = loom, dgem = dgem, chunk.size = chunk.size, display.progress = display.progress)
    # col_attrs
    print("Adding column attributes...")
    loom$create_group("col_attrs")
    add_col_attr(loom = loom, key = CA_CELLID, value = as.character(cn))
    # Check if matrix is raw counts
    if(sum(dgem%%1!=0) == 0) {
      print("Adding default metrics nUMI...")
      if(is.dgem.sparse) {
        nUMI<-Matrix::colSums(dgem)
      } else {
        nUMI<-colSums(dgem)
      }
      add_col_attr(loom = loom, key = "nUMI", value = nUMI, as.metric = T)
    } else {
      warning("Default metric nUMI was not added because the input matrix does not seem to be the raw counts.")
    }
    print("Adding default metrics nGene...")
    if(is.dgem.sparse) {
      nGene<-Matrix::colSums(dgem > 0)
    } else {
      nGene<-colSums(dgem > 0)
    }
    add_col_attr(loom = loom, key = "nGene", value = nGene, as.metric = T)
    if(!is.null(x = default.embedding)) {
      print("Adding default embedding...")
      # Add the default embedding
      if(is.null(x = default.embedding.name)) {
        default.embedding.name<-"Embedding"
      }
      add_embedding(loom = loom, embedding = as.data.frame(default.embedding), name = default.embedding.name, is.default = T)
    } else {
      warning("No default embedding set for the loom. You'll not be able to visualize it in SCope.")
    }
    # row_attrs
    print("Adding row attributes...")
    loom$create_group("row_attrs")
    add_row_attr(loom = loom, key = RA_GENE_NAME, value = as.character(rn))
    # Check if Flybase gene mapping is not empty
    if(!is.null(fbgn.gn.mapping.file.path)) {
      add_fbgn(loom = loom, dgem = dgem, fbgn.gn.mapping.file.path = fbgn.gn.mapping.file.path)
    }
    # col_edges
    print("Adding columns edges...")
    col.edges<-loom$create_group("col_edges")
    # row_edges
    print("Adding row edges...")
    row.edges<-loom$create_group("row_edges")
    # layers
    print("Adding layers...")
    layers<-loom$create_group("layers")
  }, error = function(e) {
    flush(loom = loom)
    loom$close()
    stop(e)

  }, finally = {
    flush(loom = loom)
  })
  flush(loom = loom)
}


#'@title lookup_loom
#'@description List the content of the given .loom file handler.
#'@param loom The loom file handler.
#'@export
lookup_loom<-function(loom) {
  loom$ls(recursive=TRUE)
}

#'@title open_loom
#'@description Open loom file and return a .loom file handler.
#'@param  file.path The file path to the .loom file.
#'@param mode "r+" for read/write (default). "r" for read-only. 
#'@return A loom file handler
#'@export
open_loom<-function(file.path, mode="r+") {
  return (H5File$new(file.path, mode=mode))
}



##############################
# Utils                      #
##############################

#'@title check_missings
#'@description Check if some elements are missing between given vectors x an y
#'@param x A vector of elements
#'@param y Another vector of elements
#'@return A list with "result", a logical whether the two vectors are conserved
#'                    "x.missing.in.y", a vector of elements in x that are missing in y
#'                    "y.missing.in.x", a vector of elements in y that are missing in x
check_missings<-function(x, y, check.both = TRUE) {
  x.diff.y<-setdiff(x = x, y = y)
  if(check.both) {
    y.diff.x<-setdiff(x = y, y = x)
  } else {
    y.diff.x<-NULL
  }
  return (list("result"=length(x.diff.y)==0&length(y.diff.x)==0, "x.missing.in.y"=x.diff.y,"y.missing.in.x"=y.diff.x))
}

#'@title is_order_conserved
#'@description Check if order between elements of the given vectors x an y are the same
#'@param x A vector of elements
#'@param y Another vector of elements
is_order_conserved<-function(x,y) {
  return (isTRUE(x = all.equal(current = x, target = y)))
}

#'@title open_loom
#'@description Open loom file and return a .loom file handler.
#'@param file.path The file path to the .loom file.
compress_gzb64<-function(c) {
  return (base64enc::base64encode(what = memCompress(from = c, type = "gzip")))
}

#'@title open_loom
#'@description Open loom file and return a .loom file handler.
#'@param file.path The file path to the .loom file.
decompress_gzb64<-function(gzb64c) {
  return (rawToChar(memDecompress(from = base64enc::base64decode(what = gzb64c), type = "gzip", asChar = F), multiple = F))
}

#'@title hdf5_utf8_encode
#'@description Open loom file and return a .loom file handler.
#'@param value The value from which the dtype should be UTF8 encoded.
#'@param dtype The dtype of the given value.
hdf5_utf8_encode<-function(value, dtype) {
  if(class(dtype)[1] == "H5T_STRING") {
    # Set to UTF-8 encoding otherwise read as byte array
    dtype<-dtype$set_cset('UTF-8')
  } else if(class(dtype)[1] == "H5T_COMPOUND") {
    dtypes<-unique(sapply(value, class))
    if(length(dtypes) > 1) {
      stop("Adding a data.frame as column requires the values to be the same type.")
    }
    if(dtypes == "character") {
      # Set all the columns to UTF-8 encoding otherwise read as byte array
      dtype<-H5T_COMPOUND$new(labels = colnames(value), dtypes = lapply(X = seq_along(colnames(value)), FUN = function(x) { return (H5T_STRING$new(size = Inf)$set_cset(cset = "UTF-8")) }))
    }
  }
  return (dtype)
}

###########################
# Utils
###########################

#' @title chunk_points
#' @description  Generate chunk points.
#' @author mojaveazure
#' @param data.size How big is the data being chunked.
#' @param chunk.size How big should each chunk be.
#' @return A matrix where each column is a chunk, row 1 is start points, row 2 is end points.
#' @export
chunk_points<-function(data.size, chunk.size) {
  return(vapply(
    X = 1L:ceiling(data.size / chunk.size),
    FUN = function(i) {
      return(c(
        start = (chunk.size * (i - 1L)) + 1L,
        end = min(chunk.size * i, data.size)
      ))
    },
    FUN.VALUE = numeric(length = 2L)
  ))
}


#'@title get_dtype
#'@description Get HDF5 data types.
#'@author mojaveazure
#'@param x An R object or string describing HDF5 datatype.
#'@return The corresponding HDF5 data type.
#'@import hdf5r
#'@seealso \link{hdf5r::h5types}
#'@export
get_dtype<-function(x) {
  return(switch(
    EXPR = class(x = x),
    'numeric' = h5types$double,
    'integer' = h5types$int,
    'character' = H5T_STRING$new(size = Inf),
    'logical' = H5T_LOGICAL$new(),
    stop(paste("Unknown data type:", class(x = x)))
  ))
}

#'@title get_dspace
#'@description Get the HDF5 dataspace interface object given x.
#'@param x The dataspace interface name. Either scalar or simple.
#'@return The corresponding HDF5 dataspace interface object.
get_dspace<-function(x) {
  dspaces<-c("scalar", "simple")
  if(!("scalar" %in% dspaces)) {
    stop(paste("Wrong dspace. Choose either scalar or simple."))
  }
  return (H5S$new(type = x, dims = NULL, maxdims = NULL))
}

#'@title is_base64_encoded
#'@description Checks whether the given value is base64 encoded.
#'@param value The value to check.
#'@return T is base64 encoded otherwise F.
is_base64_encoded<-function(value) {
  return (grepl(pattern = BASE64_PATTERN, x = value, perl=TRUE))
}

#'@title is_gzipped
#'@description Checks whether the given value is gzip compressed.
#'@param value The value to check.
#'@return T if gzip compressed otherwise F.
is_gzipped<-function(value) {
  tryCatch({
    memDecompress(from = value, type = "gzip", asChar = F)
    return (T)
  }, error = function(e) {
    return (F)
  }, finally = {
  })
}

#'@title is_json
#'@description Checks whether the given value is a json character
#'@param value The value to check.
#'@return T if json string character otherwise F.
is_json <-function(value) {
  tryCatch({
    rjson::fromJSON(json_str = value)
    return (T)
  }, error = function(e) {
    return (F)
  }, finally = {
  })
}

#*******************#
#*******************#
#  DATA EXTRACTION  #
#*******************#
#*******************#

#'@title get_cell_ids
#'@description Get the cell names
#'@param loom           The loom file handler.
#'@export
get_cell_ids<-function(loom
                       , is.flybase.gn = F) {
  ra<-loom[["col_attrs"]]
  return (ra[[CA_CELLID]][])
}

#'@title get_genes
#'@description Get the gene names either symbols or Flybase gene identifiers.
#'@param loom           The loom file handler.
#'@param is.flybase.gn  Whether to retrieve the Flybase gene identifiers or not.
#'@export
get_genes<-function(loom
                    , is.flybase.gn = F) {
  ra<-loom[["row_attrs"]]
  if(is.flybase.gn) {
    return (ra[["FBgn"]])
  }
  return (ra[[RA_GENE_NAME]][])
}
  
#'@title get_dgem
#'@description Get the expression matrix for the given .loom.
#'@param loom The loom file handler.
#'@return The gene expression matrix (genes as rows, cells as columns).
#'@export
get_dgem<-function(loom) {
  dgem<-t(loom[["matrix"]][,])
  cell.ids<-get_cell_ids(loom = loom)
  # Set the Cell IDs
  colnames(dgem)<-cell.ids
  genes<-get_genes(loom = loom)
  row.names(dgem)<-genes
  return (dgem)
}

#'@title get_default_embedding
#'@description Get the default embedding for the given .loom.
#'@param loom The loom file handler.
#'@return The default embedding.
#'@export
get_default_embedding<-function(loom) {
  return (loom[["col_attrs"]][[CA_EMBEDDING_NAME]])
}

#' @title get_embeddings
#' @description Get embeddings (e.g. t-SNE, UMAP) from the given loom file
#' @param loom .loom file name
#' @return List with the embeddings (each of them as data.frame)
#' @examples
#' loomPath <-  file.path(system.file('extdata', package='SCopeLoomR'), "example.loom")
#' loom <- open_loom(loomPath)
#' 
#' embeddings <- get_embeddings(loom)
#' names(embeddings)
#' head(embeddings[[1]])
#' @export
get_embeddings <- function(loom){
  embeddings <- list()
  for(embedding in get_global_meta_data(loom)$embeddings)
  {
    if(embedding$id == -1){
      tmp <- as.matrix(get_default_embedding(loom)[])
      rownames(tmp) <- get_cell_ids(loom)
      embeddings[[embedding$name]] <- tmp
    }else{
      tmp <- cbind(loom[["col_attrs"]][["Embeddings_X"]][][,as.character(embedding$id)],
                   loom[["col_attrs"]][["Embeddings_Y"]][][,as.character(embedding$id)])
      rownames(tmp) <- get_cell_ids(loom)
      embeddings[[embedding$name]] <- tmp
    }
  }
  return(embeddings)
}

#'@title get_clustering_idx_by_cluster_name
#'@description Get index of the clutering related to the given cluster.name.
#'@param loom The loom file handler.
#'@param cluster.name The name of the cluster.
#'@return The index of the clustering in the clusterings metadata global attribute corresponding to the clustering where the given cluster.name is found.
get_clustering_idx_by_cluster_name<-function(loom
                                             , cluster.name) {
  library(rlist)
  # Get global meta data
  md<-get_global_meta_data(loom = loom)
  # Unlist the nested meta data tree
  tmp<-list.flatten(md$clusterings)
  # Look for the given cluster.name
  idx = match(x = cluster.name, table = tmp)
  if(is.na(idx)) {
    return (idx)
  }
  # Reverse search from idx till the first key id is found
  for(i in idx:1) {
    if(names(tmp)[i] == "id") id<-tmp[[i]]; break
  }
  return (list.findi(.data = md$clusterings, id == id))
}

#'@title get_cluster_info_by_cluster_name
#'@description Get cluster information (Clustering ID, Clustering Name, Clustering Group, Cluster ID, Cluster Name) of the given cluster.name.
#'@param loom The loom file handler.
#'@param cluster.name The name of the cluster.
#'@return The index of the clustering in the clusterings metadata global attribute corresponding to the clustering where the given cluster.name is found.
#'@export
get_cluster_info_by_cluster_name<-function(loom
                                           , cluster.name) {
  library(rlist)
  # Get global meta data
  md<-get_global_meta_data(loom = loom)
  # Get the index of the clustering in the meta data clusterings
  clustering.idx<-get_clustering_idx_by_cluster_name(loom = loom, cluster.name = cluster.name)
  if(is.na(x = clustering.idx)) {
    stop(paste0("The given cluster ", cluster.name, " does not exist in this .loom."))
  }
  clustering.id<-md$clusterings[[clustering.idx]]$id
  cluster.idx<-list.findi(.data = md$clusterings[[clustering.idx]]$clusters, description == cluster.name)
  return (list(clustering.id=clustering.id
               , clustering.name=md$clusterings[[clustering.idx]]$name
               , clustering.group=md$clusterings[[clustering.idx]]$group
               , cluster.id=md$clusterings[[clustering.idx]]$clusters[[cluster.idx]]$id
               , cluster.name=md$clusterings[[clustering.idx]]$clusters[[cluster.idx]]$description))
}

#'@title get_clusterings
#'@description Get clustering of the given loom.
#'@param loom The loom file handler.
#'@param cluster.name The name of the cluster.
#'@return A N-by-M data.frame containing the clusterings of the given loom. N represents the cells and M the clusterings.
#'@export
get_clusterings<-function(loom) {
  return (get_col_attr_by_key(loom = loom, key = CA_CLUSTERINGS_NAME))
}

#'@title get_clustering_by_id
#'@description Get clustering with the given clustering.id of the given loom.
#'@param loom The loom file handler.
#'@param clustering.id The ID of the clustering.
#'@return A N-by-1 vector containing the cell assignments to each of the clusters of the clustering.
#'@export
get_clustering_by_id<-function(loom
                               , clustering.id) {
  ca.clusterings<-get_clusterings(loom = loom)
  return (ca.clusterings[, colnames(ca.clusterings)%in%clustering.id])
}

#' @title get_clusterings_withName
#' @description Get cell clusterings stored in the given loom file 
#' @param loom .loom file name
#' @return data.frame containing the clusterings in the given loom (cells as rows, clustering alternatives as columns).
#' @examples
#' # clusters <- get_clusterings_withName(loom)
#' # head(clusters)
#' @export
get_clusterings_withName <- function(loom){
  clusterings <- get_clusterings(loom)
  rownames(clusterings) <- get_cell_ids(loom)
  
  for(i in seq_along(colnames(clusterings))){
    clLevels <- sapply(get_global_meta_data(loom)$clusterings[[i]]$clusters, function(x) setNames(x$description, x$id))
    clusterings[,i] <- factor(clusterings[,i], levels=names(clLevels), labels=clLevels)
  }
  return(clusterings)
}

#' @title get_cellAnnotation
#' @description Get cell annotations from the given loom file
#' @param loom .loom file name
#' @return data.frame containing the cell annotations in the given loom (cells as rows, annotation variables as columns).
#' @examples
#' loomPath <-  file.path(system.file('extdata', package='SCopeLoomR'), "example.loom")
#' loom <- open_loom(loomPath)
#' 
#' cellAnnots <- get_cellAnnotation(loom)
#' head(cellAnnots)
#' @export
get_cellAnnotation <- function(loom, annotCols=NULL)  
{
  excludeCols <- c("CellID", "ClusterID", "Clusterings", "Embedding", "Embeddings_X", "Embeddings_Y", "RegulonsAUC")
  if(is.null(annotCols)) annotCols <- names(loom[["col_attrs"]])[which(!names(loom[["col_attrs"]]) %in% excludeCols)]
  
  cellAnnot <- list()
  for (colAttr in annotCols){
    cellAnnot[[colAttr]] <- get_col_attr_by_key(loom, colAttr) # loom[["col_attrs"]][[colAttr]][]
  }
  cellAnnot <- data.frame(cellAnnot)
  rownames(cellAnnot) <- get_cell_ids(loom)
  return(cellAnnot)
}


#' @title get_regulonsAuc
#' @description Get AUCell matrix from the given loom file
#' @param loom .loom file name
#' @param rows Type of data stored as rows (only for informative purposes) Default: "regulons"
#' @param columns Type of data stored as columns (only for informative purposes) Default: "cells"
#' @return AUC matrix in the slot \code{loom[["col_attrs"]][["RegulonsAUC"]][]}. 
#' If AUCell is installed, it will be stored as aucellResults object.
#' @examples
#' loomPath <-  file.path(system.file('extdata', package='SCopeLoomR'), "example.loom")
#' loom <- open_loom(loomPath)
#' 
#' regulonsAUC <- get_regulonsAuc(loom)
#' @export
get_regulonsAuc <- function(loom, rows="regulons", columns="cells") {
  mtx <- loom[["col_attrs"]][["RegulonsAUC"]][]
  rownames(mtx) <- get_cell_ids(loom)
  mtx <- t(mtx)
  names(dimnames(mtx)) <- c(rows, columns)
  
  if("AUCell" %in% rownames(installed.packages())){
    require(AUCell)
    mtx <- new("aucellResults", SummarizedExperiment::SummarizedExperiment(assays=list(AUC=mtx)))
  }
  return(mtx)
}


#' @title get_regulons
#' @description Get regulons from the given loom file
#' @param loom .loom file name
#' @param tfAsName Whether to return only the TF name (calls \code{SCENIC::getTF()}), or the regulon name as stored in the loom file.
#' @param tfSep Character used as separator for the TF name and extra values stored in the regulon name. To be passed to SCENIC::getTF()
#' @return The regulons as incidence matrix or list.
#' @examples
#' loomPath <-  file.path(system.file('extdata', package='SCopeLoomR'), "example.loom")
#' loom <- open_loom(loomPath)
#' 
#' regulons <- get_regulons(loom)
#' head(regulons)
#' @export
get_regulons <- function(loom, tfAsName=TRUE, tfSep="_"){
  incidMat <- loom[["row_attrs"]][["Regulons"]][] # incid mat
  rownames(incidMat) <- get_genes(loom)
  incidMat <- t(incidMat)
  
  if(tfAsName) rownames(incidMat) <- SCENIC::getTF(rownames(incidMat), tfSep)
  
  return(incidMat)
}

#' @title get_regulonThresholds
#' @description Get the AUC thresholds for the regulons in the given loom file
#' @param loom .loom file name
#' @param onlySelected Whether to return only the selected threshold for each regulon (as vector), 
#' or also the alternative thresholds calculated by AUCell (as list).
#' @examples
#' loomPath <-  file.path(system.file('extdata', package='SCopeLoomR'), "example.loom")
#' loom <- open_loom(loomPath)
#' 
#' thresholds <- get_regulonThresholds(loom)
#' head(thresholds)
#' @export
get_regulonThresholds <- function(loom, onlySelected=TRUE)
{
  if(onlySelected) {
    thresholds <- sapply(get_global_meta_data(loom)$regulonThresholds, function(x) setNames(x$defaultThresholdValue,x$regulon))
  }else{
    regulonNames <- sapply(get_global_meta_data(loom)$regulonThresholds, function(x) x$regulon)
    thresholds <- sapply(get_global_meta_data(loom)$regulonThresholds, function(x) {
      #regulonName <- x$regulon
      list(aucThr=list(
        selected=x$defaultThresholdValue,
        thresholds=data.frame(x$allThresholds) # if pySCENIC: transpose?
      ))
    })
  }
  
  return(thresholds)
}

#'@title get_cell_mask_by_cluster_name
#'@description Get a cell mask for the given cluster.name of the given loom.
#'@param loom The loom file handler.
#'@param cluster.name The name of the cluster.
#'@return A N-by-1 boolean vector specifying which cells belong to the given cluster.name in the given loom.
get_cell_mask_by_cluster_name<-function(loom
                                        , cluster.name) {
  # Get the cluster info given the cluster.name
  cluster.info<-get_cluster_info_by_cluster_name(loom = loom, cluster.name = cluster.name)
  # Get the clustering related to the given cluster.name
  ca.clustering<-get_clustering_by_id(loom = loom, clustering.id = cluster.info$clustering.id)
  # Create the mask
  return (ca.clustering%in%cluster.info$cluster.id)
}

#'@title get_cluster_dgem_by_name
#'@description Get a subset of the digital gene expression matrix containing only the cells in the cluster annotated by the given cluster.name.
#'@param loom The loom file handler.
#'@param cluster.name The name/description of the cluster.
#'@return A N-by-M matrix containing the gene expression levels of the cells in the cluster annotated by the given cluster.name. N represents the genes and M the cells.
#'@export
get_cluster_dgem_by_name<-function(loom
                                   , cluster.name) {
  # Get the cell mask for the given cluster.name
  mask<-get_cell_mask_by_cluster_name(loom = loom, cluster.name = cluster.name)
  dgem<-get_dgem(loom = loom)
  return (dgem[, mask])
}