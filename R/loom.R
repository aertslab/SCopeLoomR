#
# loom.R
#
# library(hdf5r)

##############################
# Global Meta data functions #
##############################

#'@title add_global_md_annotation
#'@description  Add the embedding with the given name to the global MetaData attribute.
#'@param loom   The loom file handler.
#'@param name   The name of the annotation.
#'@param values The values of the annotation to be added.
#'@export
add_global_md_annotation<-function(loom
                                 , name
                                 , values) {
  gmd<-get_global_meta_data(loom = loom)
  a<-gmd[["annotations"]]
  a[[length(a)+1]]<-list(name = name, values = as.character(unique(values)))
  gmd[["annotations"]]<-NULL
  gmd[["annotations"]]<-a
  update_global_meta_data(loom = loom, meta.data.json = rjson::toJSON(x = gmd))
}

#'@title sub_embedding_group_exists
#'@description  Check whether a sub embedding group of the given parent.clustering.id exists in the MetaData
#'@param parent.clustering.id The id of the parent clustering
#'@export
sub_embedding_group_exists<-function(parent.clustering.id) {
  gmd<-get_global_meta_data(loom = loom)
  e<-gmd[["embeddings"]]
  x<-sapply(e, function(x) { x['parentClusteringId'] == parent.clustering.id })
  return (is.na(match(TRUE,x)))
}


#'@title init_global_md_embeddings
#'@description  Remove all embeddings from the given loom
#'@param loom   The loom file handler.
#'@export
init_global_md_embeddings<-function(loom) {
  gmd<-get_global_meta_data(loom = loom)
  gmd[["embeddings"]]<-NULL
  gmd[["embeddings"]]<-list()
  update_global_meta_data(loom = loom, meta.data.json = rjson::toJSON(x = gmd))
}

#'@title add_global_md_embedding
#'@description  Add the embedding with the given name to the global MetaData attribute.
#'@param loom       The loom file handler.
#'@param name       The name of the embedding to add.
#'@param is.default Is the given embedding the default embedding to use in the .loom file.
#'@export
add_global_md_embedding<-function(loom
                                , id
                                , name
                                , is.default = F
                                , parent.clustering.id = NULL
                                , cluster.id = NULL) {
  gmd<-get_global_meta_data(loom = loom)
  e<-gmd[["embeddings"]]
  if(is.default) {
    default<-"true"
  } else {
    default<-"false"
  }
  
  if(is.null(parent.clustering.id)) {
    parent.clustering.id<-(-1) # By default set to all cells
  }
  if(is.null(cluster.id)) {
    cluster.id<-(-1) # By default set to all cells
  }
  
  e[[length(e)+1]]<-list(id = id, name = name, parentClusteringId = parent.clustering.id, clusterId = cluster.id, default = default)
  gmd[["embeddings"]]<-NULL
  gmd[["embeddings"]]<-e
  update_global_meta_data(loom = loom, meta.data.json = rjson::toJSON(x = gmd))
}


#'@title sub_clustering_exists
#'@description  Check whether a sub clustering of the given parent.clustering.id exists in the MetaData
#'@param parent.clustering.id The clustering id of the parent
#'@export
sub_clustering_exists<-function(parent.clustering.id) {
  gmd<-get_global_meta_data(loom = loom)
  c<-gmd[["clusterings"]]
  x<-sapply(c, function(x) { x['parentClusteringId'] == parent.clustering.id })
  return (is.na(match(TRUE,x)))
}

#'@title add_global_md_clustering
#'@description  Add the clustering annotation to the global MetaData attribute.
#'@param loom                               The loom file handler.
#'@param group                              The name of the group of clusterings
#'@param name                               The name given to the given clustering.
#'@param clusters                           A list of the the cluster id for each cell present in the same order as in the columns of gene expression matrix
#'@param annotation.cluster.id.cl           The column name to use for the IDs of the clusters found by the given clustering group.
#'@param annotation.cluster.description.cl  The column name to use for the description of the clusters found by the given clustering group.
#'@export
add_global_md_clustering<-function(loom
                                 , id 
                                 , level = 0
                                 , group
                                 , name
                                 , params
                                 , parent.clustering.id = NULL
                                 , clusters
                                 , annotation = NULL
                                 , annotation.cluster.id.cl = NULL
                                 , annotation.cluster.description.cl = NULL) {
  gmd<-get_global_meta_data(loom = loom)
  c<-gmd[["clusterings"]]
  unique.clusters<-sort(as.numeric(unique(clusters))-1, decreasing = F)
  clusters<-lapply(X = unique.clusters, FUN = function(cluster.id) {
    description<-paste("NDA - Cluster", cluster.id)
    if(!is.null(annotation)) {
      # If annotation for the current cluster not empty then add
      d<-annotation[annotation[[annotation.cluster.id.cl]] == cluster.id, annotation.cluster.description.cl]
      if(nchar(d)>0) {
        description<-paste0(d, " (",cluster.id,")")
      }
    }
    return (list(id = cluster.id
                 , description = description))
  })
  ca<-loom[["col_attrs"]]
  clusterings<-ca[["Clusterings"]][]
  # Check if sub-clustering
  if(is.null(parent.clustering.id)) {
    parent.clustering.id<-(-1)
  }
  clustering<-list(id = id
                 , level = level
                 , parentClusteringId = parent.clustering.id
                 , group = group
                 , name = name
                 , params = params
                 , clusters = clusters)
  c[[length(c)+1]]<-clustering
  gmd[["clusterings"]]<-NULL
  gmd[["clusterings"]]<-c
  update_global_meta_data(loom = loom, meta.data.json = rjson::toJSON(x = gmd))
  loom$flush
}

#'@title add_global_md_regulon_thresholds
#'@description  Add the regulon thresholds annotation to the global MetaData attribute.
#'@param loom                           The loom file handler.
#'@param regulon.threshold.assignments  The automated regulon cell assignments object generated by SCENIC.
#'@param regulon.enrichment.table       The regulon enrichment table generated by SCENIC.
#'@export
add_global_md_regulon_thresholds<-function(loom
                                         , regulon.threshold.assignments
                                         , regulon.enrichment.table) {
  gmd<-get_global_meta_data(loom = loom)
  rT<-gmd[["regulonThresholds"]]
  for(regulon in names(regulon.threshold.assignments)) {
    regulon.name<-strsplit(regulon, "_")[[1]][1]
    # Store the threshold assignments
    rta<-regulon.threshold.assignments[[regulon]]
    AUC.thresholds<-rta$aucThr$thresholds
    allThresholds<-do.call(what = "c", args = lapply(X = row.names(AUC.thresholds), FUN = function(threshold.name) {
      l<-list()
      l[[threshold.name]]<-AUC.thresholds[threshold.name,"threshold"]
      return (l)
    }))
    AUC.selected.threshold<-rta$aucThr$selected
    # Store the motif name
    if(!("NES"%in%colnames(regulon.enrichment.table))) {
      stop("Error: NES column is not provided in the regulon enrichment table!")
    }
    if(!("motif"%in%colnames(regulon.enrichment.table))) {
      stop("Error: motif column is not provided in the regulon enrichment table!")
    }
    regulon.enrichment.table<-regulon.enrichment.table[with(regulon.enrichment.table, order(NES, decreasing = T)), ]
    motif.name = regulon.enrichment.table$motif[regulon.enrichment.table$highlightedTFs == regulon.name][1]
    regulon.tresholds<-list(regulon = gsub(pattern = " ", replacement = "_", x = regulon)
                            , defaultThresholdValue = as.numeric(AUC.selected.threshold)
                            , defaultThresholdName = names(AUC.selected.threshold)
                            , allThresholds = allThresholds
                            , motifData = paste0(motif.name,".png"))
    rT[[length(rT)+1]]<-regulon.tresholds
  }
  gmd[["regulonThresholds"]]<-NULL
  gmd[["regulonThresholds"]]<-rT
  loom$flush
  update_global_meta_data(loom = loom, meta.data.json = rjson::toJSON(x = gmd))
  loom$flush
}

#'@title get_global_meta_data
#'@description Get the global MetaData attribute as a R object.
#'@param loom The loom file handler.
#'@export
get_global_meta_data<-function(loom) {
  meta.data<-decompress_gzb64(gzb64c = h5attr(x = loom, which = "MetaData"))
  return (rjson::fromJSON(json_str = meta.data))
}

update_global_meta_data<-function(loom
                                  , meta.data.json) {
  compressed.meta.data<-compress_gzb64(c = as.character(meta.data.json))
  h5attr(x = loom, which = "MetaData")<-compressed.meta.data
}

#'@export
init_global_meta_data<-function(loom) {
  meta.data<-list(annotations = list()
                , embeddings = list()
                , clusterings = list()
                , regulonThresholds = list())
  meta.data.json<-rjson::toJSON(meta.data)
  if(!("MetaData"%in%list.attributes(object = loom))) {
    compressed.meta.data<-compress_gzb64(c = as.character(meta.data.json))
    loom$create_attr(attr_name = "MetaData", robj = as.character(compressed.meta.data), dtype = H5T_STRING$new(size=Inf))
  } else {
    update_global_meta_data(loom = loom, meta.data.json = as.character(meta.data.json))
  }
}

#######################
# Embedding functions #
#######################

sub_embedding_exists<-function(loom, id, embedding) {
  # Only test on X (it should be true for Y)
  ca.embeddings.x<-get_col_attr_by_key(loom = loom, key = "Embeddings_X")
  # Retrieve sub embedding groups with the given id
  ca.embeddings.x.ss<-ca.embeddings.x[, grepl(pattern = paste0("^",id,"_[0-9]+"), x = colnames(ca.embeddings.x))]
  # By default set the embedding group id to store the given embedding to N (start at 0) 
  se.group.id<-ncol(ca.embeddings.x.ss)
  for(i in 1:ncol(ca.embeddings.x.ss)) {
    tmp<-ca.embeddings.x.ss[,i]
    # Check if cell embeddings are set in the current embedding group
    nb.es<-tmp[cell.ids%in%row.names(embedding)] == -1
    if(sum(nb.es) == nrow(embedding)) {
      return (i)
    }
  }
  return (se.group.id)
}

add_default_embedding<-function(loom, embedding) {
  embedding<-as.data.frame(embedding)
  colnames(embedding)<-c("_X","_Y")
  add_col_attr(loom = loom, key = "Embedding", value = embedding)
}

append_embedding_group_update_ca<-function(loom
                                      , embedding
                                      , embedding.group.id
                                      , coord.labels) {
  cell.ids<-get_cells(loom = loom)
  for(i in seq_along(coord.labels)) {
    ca.embeddings<-get_col_attr_by_key(loom = loom, key = coord.labels[i])
    sub.embedding<-data.frame(x = rep(-1, nrow(ca.embeddings)))
    sub.embedding[cell.ids%in%row.names(embedding)]<-embedding[,i]
    # Update the embedding
    colnames(e)<-embedding.group.id
    ca.embeddings<-cbind(ca.embeddings, e)
    # Update the current coordinates Embeddings
    update_col_attr(loom = loom, key = coord.labels[i], value = as.data.frame(ca.embeddings))
  }
}

#'@title add_embedding
#'@description Add the given embedding as a row attribute and meta data related to the given embeddding to the given .loom file handler.
#'@param loom       The loom file handler.
#'@param embedding  A M-by-2 data.frame of the embeddings with M cells.
#'@param name       The name of the given embedding.
#'@param is.default Default embedding to use in the .loom file.
#'@export
add_embedding<-function(loom
                        , embedding
                        , name
                        , is.default = F
                        , parent.clustering.id = NULL
                        , parent.clustering.id.lookup.query = NULL # list("level"=0, "group"="Seurat", "param.name"="", "param.value"="")
                        , cluster.id = NULL) {
  coord.labels<-c("Embeddings_X", "Embeddings_Y")
  # Check if default and not a sub-cluster embedding
  # Add the default embedding also to Embeddings_X and Embeddings_Y
  if(is.default & is.null(parent.clustering.id) & is.null(cluster.id)) {
    add_default_embedding(loom = loom, embedding = embedding)
  }
  
  ca<-loom[["col_attrs"]]
  
  if(!is.null(parent.clustering.id)) {
    # Add a sub embedding
    if(is.null(parent.clustering.id)) {
      # Get parent clustering id based on the given look-up query
      parent.clustering.id<-get_clid_by_clustering_param(loom = loom
                                                       , clustering.level = parent.clustering.id.lookup.query["level"]
                                                       , clustering.group = parent.clustering.id.lookup.query["group"]
                                                       , clustering.param.name = parent.clustering.id.lookup.query["param.name"]
                                                       , clustering.param.value = parent.clustering.id.lookup.query["param.value"])
    }
    id<-paste0(parent.clustering.id,"_sub")
    cell.ids<-get_cells(loom = loom)
    # Check if the sub embedding already exists given the parent.clustering.id
    if(sub_embedding_group_exists(parent.clustering.id)) {
      # Get the next sub embedding group id slot
      se.group.suffix<-sub_embedding_exists(loom = loom, id = id, embedding = embedding)
      id<-paste0(id,"_",se.group.suffix)
    } else {
      id<-paste0(id,"_0")
    }
    append_embedding_group_update_ca(loom = loom, embedding = embedding, embedding.group.id = id, coord.labels = coord.labels)
  } else {
    # Add a main embedding
    if(sum(coord.labels%in%names(ca)) != 2) {
      for(i in seq_along(coord.labels)) {
        e<-as.data.frame(embedding[,i])
        id<-"0"
        colnames(e)<-id
        add_col_attr(loom = loom, key = coord.labels[i], value = e)
      }
    } else {
      for(i in seq_along(coord.labels)) {
        ca.embeddings<-get_col_attr_by_key(loom = loom, key = coord.labels[i])
        e<-as.data.frame(embedding[,i])
        id<-as.character(ncol(ca.embeddings))
        colnames(e)<-id
        ca.embeddings<-cbind(ca.embeddings, e)
        # Update the current coordinates Embeddings
        update_col_attr(loom = loom, key = coord.labels[i], value = as.data.frame(ca.embeddings))
      }
    }
  }
  loom$flush()
  # Add the given embedding to the global attribute MetaData
  add_global_md_embedding(loom = loom
                          , id = id
                          , name = name
                          , is.default = is.default
                          , parent.clustering.id = parent.clustering.id
                          , cluster.id = cluster.id)
  loom$flush()
}

#########################
# Clusterings functions #
#########################

get_clid_by_clustering_param<-function(loom
                                         , clustering.level
                                         , clustering.group
                                         , clustering.param.name
                                         , clustering.param.value) {
  ca<-loom[["col_attrs"]]
  ca.clusterings<-ca[['Clusterings']][]
  gmd<-get_global_meta_data(loom = loom)
  c<-gmd[["clusterings"]]
  x<-sapply(c, function(x) { 
    x['params'][[clustering.param.name]] == clustering.param.value &
      x['level'] == clustering.level &
      x['group'] == clustering.group
  })
  clid<-colnames(ca.clusterings)[x]
  if(length(clid) != 1) {
    stop(paste0("The clustering with param '", paste(clustering.param.name, clustering.param.value),"' do not exist."))
  }
}

#'@title add_seurat_clustering
#'@description Add all the Seurat clusterings in the given seurat object to the given .loom file handler.
#'@param loom                               The loom file handler.
#'@param seurat                             The Seurat object
#'@param default.clustering.resolution      The clustering resolution (i.e.: res.2, ...) of the clustering that should be set as the default.
#'@param seurat.markers.file.path.list      The named list of file paths to the markers saved in RDS format. The names should be the resolution id of the corresponding clustering (e.g.: res2.0).
#'@param annotation                         A data.frame with annotation for the clusters
#'@param annotation.cluster.id.cl           The column name to use for the IDs of the clusters found by the given clustering group.
#'@param annotation.cluster.description.cl  The column name to use for the description of the clusters found by the given clustering group.
#'@export
add_seurat_clustering<-function(loom
                                , seurat
                                , default.clustering.resolution
                                , seurat.markers.file.path.list
                                , seurat.clustering.level
                                , parent.seurat.clustering.resolution
                                , annotation
                                , annotation.cluster.id.cl
                                , annotation.cluster.description.cl) {
  # Get parent clustering id based on parent Seurat clustering resolution
  parent.seurat.clustering.id<-get_clid_by_clustering_param(loom = loom
                                                          , clustering.level = seurat.clustering.level
                                                          , clustering.group = "Seurat"
                                                          , clustering.param.name = "resolution"
                                                          , clustering.param.value = parent.seurat.clustering.resolution)
  clustering.resolutions<-as.numeric(stringr::str_split_fixed(string = names(seurat@calc.params)[grep(pattern = "FindClusters", x = names(seurat@calc.params))], pattern = "FindClusters.res.", n = 2)[,2])
  for(res in clustering.resolutions) {
    seurat<-seurat::SetAllIdent(object=seurat, id=paste0("res.",res))
    cluster.ids<-seurat@ident
    is.default.clustering<-F
    # Add the Seurat clusters
    print("Adding Seurat clusters...")
    a<-NULL
    ac.id.cl<-NULL
    ac.description.cl<-NULL
    if(res == default.clustering.resolution) {
      print("Adding default Seurat clusters...")
      a<-annotation
      ac.id.cl<-annotation.cluster.id.cl
      ac.description.cl<-annotation.cluster.description.cl
      is.default.clustering<-T
    }
    loom$flush()
    clid<-add_clustering(loom = loom
                       , group = "Seurat"
                       , name = paste("Seurat, resolution",res)
                       , params = list("resolution"=res)
                       , clusters = cluster.ids
                       , clustering.level = seurat.clustering.level
                       , parent.clustering.id = parent.seurat.clustering.id
                       , is.default = is.default.clustering
                       , annotation = a
                       , annotation.cluster.id.cl = ac.id.cl
                       , annotation.cluster.description.cl = ac.description.cl)
    loom$flush()
    # Add the Seurat markers
    print("Adding Seurat markers...")
    seurat.markers<-readRDS(file = seurat.markers.file.path.list[[paste0("res.", res)]])
    seurat.m.list<-split(x = seurat.markers, f = seurat.markers$cluster)
    seurat.m.list<-lapply(X = seurat.m.list, function(cluster) {
      return (cluster$gene)
    })
    add_clustering_markers(loom = loom, clustering.id = clid, clustering.markers = seurat.m.list)
    loom$flush()
  }
}

append_clustering_update_ca<-function(loom
                                      , clustering.id
                                      , clustering) {
  k<-"Clusterings"
  ca.clusterings<-get_col_attr_by_key(loom = loom, key = k)
  colnames(clustering)<-clustering.id
  # Append this clustering
  ca.clusterings<-cbind(ca.clusterings, clustering)
  update_col_attr(loom = loom, key = k, value = as.data.frame(x = ca.clusterings))
}

#'@title add_clustering
#'@description Add the given clusters in the given group column attribute and meta data related to the given clustering to the given .loom file handler.
#'@param loom                               The loom file handler.
#'@param group                              The for the given clustering group to which the given clusters have to be added
#'@param name                               The name given to this clustering
#'@param parent.clustering.id               Clustering ID of the parent in case it's a subclustering.
#'@param clusters                           A named list of the cell id and assigned the cluster id.
#'@param is.default                         Set this clustering be set as default one.
#'@param annotation                         A data.frame with annotation for the clusters
#'@param annotation.cluster.id.cl           The column name to use for the IDs of the clusters found by the given clustering group.
#'@param annotation.cluster.description.cl  The column name to use for the description of the clusters found by the given clustering group.
#'@export
add_clustering<-function(loom
                         , group
                         , name
                         , params
                         , clusters
                         , level = 0
                         , parent.clustering.id = NULL
                         , is.default = F
                         , annotation = NULL
                         , annotation.cluster.id.cl =  NULL
                         , annotation.cluster.description.cl = NULL) {
  id<-0
  # If the clustering is the default one
  # Add it as the generic column attributes ClusterID and ClusterName
  if(is.default) {
    add_col_attr(loom = loom, key = "ClusterID", value = clusters)
    unique.clusters<-sort(as.numeric(unique(clusters))-1, decreasing = F)
    for(cluster in unique.clusters) {
      description<-paste0("NDA - Cluster ", cluster)
      names(clusters)<-description
      if(!is.null(annotation)) {
        description<-annotation[annotation[[annotation.cluster.id.cl]] == cluster, annotation.cluster.description.cl]
      }
      names(clusters)[clusters == cluster]<-description
    }
    add_col_attr(loom = loom, key = "ClusterName", value = names(clusters))
  }
  # Adding the clustering data
  k<-"Clusterings"
  if(col_attrs_exists_by_key(loom = loom, key = k)) {
    print("Clusterings already exists...")
    ca.clusterings<-get_col_attr_by_key(loom = loom, key = k)
    # Set the clustering id
    id<-ncol(ca.clusterings)+1 # n clusterings
    # Check if we are adding a sub clustering
    if(!is.null(parent.clustering.id)) {
      id<-paste0(parent.clustering.id,"_sub")
      cell.ids<-get_cells(loom = loom)
      # Check if the sub clustering already exists given the parent.clustering.id
      if(sub_clustering_exists(parent.clustering.id)) {
        ca.clusterings[[id]][cell.ids%in%names(clusters)]<-clusters
        update_col_attr(loom = loom, key = k, value = as.data.frame(x = ca.clusterings))
      } else {
        sub.clustering<-data.frame(x = rep(-1, nrow(ca.clusterings)))
        sub.clustering['x'][cell.ids%in%names(clusters)]<-clusters
        append_clustering_update_ca(loom = loom, clustering.id = id, clustering = sub.clustering)
      }
    } else {
      clustering<-data.frame(x = clusters)
      append_clustering_update_ca(loom = loom, clustering.id = id, clustering = clustering)
    }
  } else {
    print("Clusterings created...")
    clustering<-data.frame("0" = clusters)
    add_col_attr(loom = loom, key = k, value = as.data.frame(x = clustering))
  }
  loom$flush()
  # Adding the clustering meta data
  add_global_md_clustering(loom = loom
                           , id = id
                           , level = level
                           , group = group
                           , name = name
                           , params = params
                           , clusters = clusters
                           , parent.clustering.id = parent.clustering.id
                           , annotation = annotation
                           , annotation.cluster.id.cl = annotation.cluster.id.cl
                           , annotation.cluster.description.cl = annotation.cluster.description.cl)
  loom$flush()
  return (id)
}

####################
# SCENIC functions #
####################

#'@title add_scenic_regulons
#'@description Add the regulons with their target genes generated by SCENIC as a row attribute to the given .loom file handler.
#'@param loom     The loom file handler.
#'@param dgem     A matrix of the gene expression with M genes as rows and N cells as columns.
#'@param regulons A list of list of the regulons and their target genes generated by SCENIC.
#'@export
add_scenic_regulons<-function(loom
                              , dgem
                              , regulons) {
  regulons.mask<-do.call(what = "cbind", args = lapply(seq_along(regulons), function(regulon.idx) {
    reg.name<-names(regulons)[regulon.idx]
    reg.genes<-regulons[regulon.idx]
    reg.mask<-data.frame("x"=row.names(dgem)%in%reg.genes, stringsAsFactors = F)
    colnames(reg.mask)<-reg.name
    return (reg.mask)
  }))
  row.names(regulons.mask)<-row.names(dgem)
  colnames(regulons.mask)<-gsub(pattern = " ", replacement = "_", x = colnames(regulons.mask))
  add_row_attr(loom = loom, key = "Regulons", value = as.data.frame(x = regulons.mask))
  loom$flush()
}

#'@title add_scenic_regulons_auc_matrix
#'@description Add the regulons AUC matrix generated by SCENIC as a column attribute to the given .loom file handler.
#'@param loom         The loom file handler.
#'@param regulons.AUC A matrix of the regulons AUC values with M regulons as rows and N cells as columns.
#'@export
add_scenic_regulons_auc_matrix<-function(loom
                                         , regulons.AUC) {
  add_col_attr(loom = loom, key = "RegulonsAUC", value = as.data.frame(x = t(regulons.AUC)))
  loom$flush()
}

###########################
# Col Meta data functions #
###########################

#'@title get_cells
#'@description Get the cell names
#'@param loom           The loom file handler.
#'@export
get_cells<-function(loom
                    , is.flybase.gn = F) {
  ra<-loom[["col_attrs"]]
  return (ra[["CellID"]][])
}

###########################
# Row Meta data functions #
###########################

#'@title add_clustering_markers
#'@description Add the clustering markers as a row attribute to the given .loom file handler.
#'@param loom           The loom file handler.
#'@param dgem           A matrix of the gene expression with M genes as rows and N cells as columns.
#'@param clustering.id  The clustering id that the given clustering.markers are specific for.
#'@param regulons       A list of list of the clustering markers.
#'@export
add_clustering_markers<-function(loom
                               , clustering.id
                               , clustering.markers) {
  genes<-get_genes(loom = loom, is.flybase.gn = F)
  clustering.markers.mask<-do.call(what = "cbind", args = lapply(seq_along(clustering.markers), function(cluster.idx) {
    cluster.name<-names(clustering.markers)[cluster.idx]
    cluster.markers<-clustering.markers[cluster.idx]
    cm.mask<-data.frame("x" = genes %in% cluster.markers, stringsAsFactors = F)
    colnames(cm.mask)<-cluster.name
    return (cm.mask)
  }))
  row.names(clustering.markers.mask)<-genes
  clustering.id<-ncol(loom[["col_attrs"]][["Clusterings"]][])-1
  print(paste0("Adding markers for clustering ID ", clustering.id))
  add_row_attr(loom = loom, key = paste0("ClusteringMarkers_",clustering.id), value = as.data.frame(x = clustering.markers.mask))
  loom$flush()
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
  return (ra[["Gene"]][])
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
  fbgn.gn.mapping<-utils::read.table(file = fbgn.gn.mapping.file.path, header = F, sep = "\t", quote = '', stringsAsFactors = F)
  colnames(fbgn.gn.mapping)<-c("FBgn","Gene")
  genes<-merge(x = data.frame("Gene"=row.names(dgem)), y = fbgn.gn.mapping, by = "Gene")
  add_row_attr(loom = loom, key = "FBgn", value = genes$FBgn)
}

#####################
# Generic functions #
#####################
#'@export
lookup_all_global_attr<-function(loom) {
  list.attributes(object = loom)
}


#'@export
remove_global_attr<-function(loom
                             , key) {
  loom$attr_delete(attr_name = key)
  loom$flush()
}


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
#'@export
add_global_attr<-function(loom
                          , key
                          , value) {
  loom$create_attr(attr_name = key, robj = value, dtype = getDtype(x = value))
  loom$flush()
}

#'@export
remove_row_attr<-function(loom
                          , key) {
  loom$link_delete(name = paste0("row_attrs/", key))
  loom$flush()
}

#'@export
update_row_attr<-function(loom
                          , key
                          , value) {
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
                       , value) {
  ra<-loom[["row_attrs"]]
  ra[[key]]<-value
  loom$flush()
}

#'@export
col_attrs_exists_by_key<-function(loom
                                , key) {
  ca<-loom[["col_attrs"]]
  return (key %in% names(ca))
}

#'@export
get_col_attr_by_key<-function(loom
                            , key) {
  ca<-loom[["col_attrs"]]
  return (ca[[key]][])
}

#'@export
remove_col_attr<-function(loom
                          , key) {
  loom$link_delete(name = paste0("col_attrs/", key))
  loom$flush()
}

#'@export
update_col_attr<-function(loom
                          , key
                          , value) {
  loom$link_delete(name = paste0("col_attrs/",key))
  add_col_attr(loom = loom, key = key, value = value)
}

#'@title add_col_attr
#'@description Add a new column attribute to the given .loom object accessible by the given key and containing the given value.
#'@param loom   The loom file handler.
#'@param key    The name of the new added attribute.
#'@param value  The value of the new added attribute.
#'@param as.md.annotation Whether to show this attribute in the compare datasets tab
#'@export
add_col_attr<-function(loom
                       , key
                       , value
                       , as.md.annotation = F) {
  ca<-loom[["col_attrs"]]
  ca[[key]]<-value
  loom$flush()
  if(as.md.annotation) {
    add_global_md_annotation(loom = loom, name = key, values = value)
    loom$flush()
  }
}

#'@title add_matrix
#'@description Add the given gene expression matrix dgem to the given .loom object.
#'@param loom The loom file handler.
#'@param dgem A matrix of the gene expression with M genes as rows and N cells as columns.
#'@export
add_matrix<-function(loom
                     , dgem
                     , chunk.size
                     , display.progress) {
  row.names(dgem)<-NULL
  colnames(dgem)<-NULL
  dtype<-getDtype(x = dgem[1, 1])
  loom$create_dataset(
    name = 'matrix',
    dtype = dtype,
    dims = rev(x = dim(x = dgem))
  )
  chunk.points<-chunkPoints(
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
  loom$flush()
}

#'@title finalize
#'@description finalize
#'@loom The loom file handler.
#'@export
finalize<-function(loom) {
  loom$flush()
  loom$close_all()
}

#'@title build_loom
#'@description build_loom
#'@param file.name                  A string naming the .loom file to be generated.
#'@param dgem                       A matrix of the gene expression with M genes as rows and N cells as columns.
#'@param default.embedding          A M-by-2 data.frame of the embedding (X and Y coordinates) of the cells.
#'@param fbgn.gn.mapping.file.path  A N-by-2 data.frame containing the mapping between the Flybase gene and the gene symbol.
#'@export
build_loom<-function(file.name
                     , title = NULL
                     , genome = NULL
                     , dgem
                     , default.embedding
                     , default.embedding.name
                     , fbgn.gn.mapping.file.path = NULL
                     , chunk.size = 1000
                     , display.progress = T) {
  loom<-H5File$new(filename = file.name, mode = "w")
  tryCatch({
    # title
    if(!is.null(title)) {
      add_global_attr(loom = loom, key = "title", value = as.character(title))
    }
    # Genome
    if(!is.null(genome)) {
      add_global_attr(loom = loom, key = "Genome", value = as.character(genome))
    }
    cn<-colnames(dgem)
    rn<-row.names(dgem)
    print("Adding global attributes...")
    # global MetaData attribute
    init_global_meta_data(loom = loom)
    # matrix
    # Check the type of the sparse matrix
    # convert to dgCMatrix if necessary to speedup populating the matrix slot
    if(class(dgem) == "dgTMatrix") {
      print("Converting to dgCMatrix...")
      dgem<-methods::as(object = dgem, Class = "dgCMatrix")
    }
    print("Adding matrix...")
    add_matrix(loom = loom, dgem = dgem, chunk.size = chunk.size, display.progress = display.progress)
    # col_attrs
    print("Adding column attributes...")
    loom$create_group("col_attrs")
    add_col_attr(loom = loom, key = "CellID", value = as.character(cn))
    print("Adding default embedding...")
    # Add the default embedding
    add_embedding(loom = loom, embedding = as.data.frame(default.embedding), name = default.embedding.name, is.default = T)
    # row_attrs
    print("Adding row attributes...")
    loom$create_group("row_attrs")
    add_row_attr(loom = loom, key = "Gene", value = as.character(rn))
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
    loom$flush()
    loom$close()
    stop(e)

  }, finally = {
    loom$flush()
  })
  loom$flush()
}

#'@export
lookup_loom<-function(loom) {
  loom$ls(recursive=TRUE)
}

#'@export
open_loom<-function(file.path) {
  return (H5File$new(file.path, mode="r+"))
}

##############################
# Utils                      #
##############################

compress_gzb64<-function(c) {
  return (base64enc::base64encode(what = memCompress(from = c, type = "gzip")))
}

decompress_gzb64<-function(gzb64c) {
  return (rawToChar(memDecompress(from = base64enc::base64decode(what = gzb64c), type = "gzip", asChar = F), multiple = F))
}

###########################
# Utils (loomR)
###########################

#' @title chunkPoints
#' @description  Generate chunk points
#' @param data.size How big is the data being chunked
#' @param chunk.size How big should each chunk be
#' @return A matrix where each column is a chunk, row 1 is start points, row 2 is end points
#' @export
chunkPoints<-function(data.size, chunk.size) {
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


#'@title getDtype
#'@description Get HDF5 data types
#'@param x An R object or string describing HDF5 datatype
#'@return The corresponding HDF5 data type
#'@ rdname getDtype
#'@import hdf5r
#'@seealso \link{hdf5r::h5types}
#'@export
getDtype<-function(x) {
  return(switch(
    EXPR = class(x = x),
    'numeric' = h5types$double,
    'integer' = h5types$int,
    'character' = H5T_STRING$new(size = Inf),
    'logical' = H5T_LOGICAL$new(),
    stop(paste("Unknown data type:", class(x = x)))
  ))
}
