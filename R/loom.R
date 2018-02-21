# loom.R
#

##############################
# Global Meta data functions #
##############################

#'
#'@description  Add the embedding with the given name to the global MetaData attribute.
#'@param loom   The loom file handler.
#'@param name   The name of the annotation.
#'@param values The uniques values of the annotation to be added.
#'
add_global_md_annotation<-function(loom
                                 , name
                                 , values) {
  gmd<-get_global_meta_data(loom = loom)
  a<-gmd["annotations"]
  append(x = a, values = list(name = name, values = unique(values)))
  gmd["annotations"]<-a
  update_global_meta_data(loom = loom, meta.data.json = rjson::toJSON(x = gmd))
  loom$flush
}

#'
#'@description  Add the embedding with the given name to the global MetaData attribute.
#'@param loom       The loom file handler.
#'@param name       The name of the embedding to add.
#'@param is.default Is the given embedding the default embedding to use in the .loom file.
#'
add_global_md_embedding<-function(loom
                                , name
                                , is.default = F) {
  gmd<-get_global_meta_data(loom = loom)
  e<-gmd["embeddings"]
  if(is.default) {
    id<-(-1)
  } else {
    id<-len(e)+1
  }
  append(x = e, values = list(id = id, name = name))
  gmd["embeddings"]<-e
  update_global_meta_data(loom = loom, meta.data.json = rjson::toJSON(x = gmd))
  loom$flush
}

#'
#'@description  Add the clustering annotation to the global MetaData attribute.
#'@param loom                               The loom file handler.
#'@param group                              The name of the group of clusterings
#'@param name                               The name given to the givem clustering.
#'@param annotation.cluster.id.cl           The column name to use for the IDs of the clusters found by the given clustering method.
#'@param annotation.cluster.description.cl  The column name to use for the description of the clusters found by the given clustering method.
#'
add_global_md_clustering<-function(loom
                                 , group
                                 , name
                                 , annotation
                                 , annotation.cluster.id.cl # Column label containing the id of the cluster
                                 , annotation.cluster.description.cl) { # Column label containing the description of the cluster
  gmd<-get_global_meta_data(loom = loom)
  c<-gmd["clusterings"]
  clusters<-do.call("c", apply(X = annotation, MARGIN = 1, FUN = function(cluster) {
    description<-cluster[[annotation.cluster.description.cl]]
    # If no annotation provided for the current cluster then just concat "Cluster" with the cluster ID
    if(nchar(cluster[[annotation.cluster.description.cl]])>0) {
      description<-paste("Cluster",cluster[[annotation.cluster.id.cl]])
    }
    return (list(id = cluster[[annotation.cluster.id.cl]]
                 , description = cluster[[annotation.cluster.description.cl]]))
  }))
  ca<-file.h5[["col_attrs"]]
  clusterings<-ca[[paste0(method,"Clusterings")]][]
  clustering<-list(id = ncol(clusterings)+1, # n clusterings
                   group = group,
                   name = name,
                   clusters = clusters)
  append(x = c, values = clustering)
  gmd["clusterings"]<-e
  update_global_meta_data(loom = loom, meta.data.json = rjson::toJSON(x = gmd))
  loom$flush
}

#'
#'@description  Add the markers to the global MetaData attribute.
#'@param loom             The loom file handler.
#'@param clustering.group The name of the group of clusterings.
#'@param clustering.id    The id of the clustering from which the markers have been calculated.
#'@param cluster.id       The id of the cluster to which the markers are specific for
#'
add_global_md_cluster_markers<-function(loom
                                      , clustering.group
                                      , clustering.id
                                      , cluster.id) {
  gmd<-get_global_meta_data(loom = loom)
  m<-gmd["markers"]
  append(x = m, values = list(clusteringGroup = clustering.group
                            , clusteringId = clustering.id
                            , clusterId = cluster.id))
  gmd["markers"]<-m
  update_global_meta_data(loom = loom, meta.data.json = rjson::toJSON(x = gmd))
  loom$flush
}

#'
#'@description Get the global MetaData attribute as a R object.
#'@param loom The loom file handler.
#'
get_global_meta_data<-function(loom) {
  return (rjson::fromJSON(json_str = h5attr(x = loom, which = "MetaData")))
}

update_global_meta_data<-function(loom
                                  , meta.data.json) {
  h5attr(x = loom, which = "MetaData")<-meta.data.json
}

add_global_meta_data<-function(loom) {
  meta.data<-list(annotations = list()
                , embeddings = list()
                , clusterings = list()
                , markers = list())
  meta.data.json<-rjson::toJSON(meta.data)
  if(!("MetaData"%in%list.attributes(object = loom))) {
    loom$create_attr(attr_name = "MetaData", robj = meta.data.json)
  } else {
    update_global_meta_data(loom = loom, meta.data.json = meta.data.json)
  }
}

#####################
# General functions #
#####################

#'
#'@description Add the given embedding as a row attribute and meta data related to the given embeddding to the given .loom file handler.
#'@param loom       The loom file handler.
#'@param embedding  A M-by-2 data.frame of the embeddings with M cells.
#'@param name       The name of the given embedding.
#'@param is.default Is the given embedding the default embedding to use in the .loom file.
#'
add_embedding<-function(loom
                        , embedding
                        , name
                        , is.default = F) {
  coord.labels<-c("Embeddings_X", "Embeddings_Y")
  if(is.default) {
    add_col_attr(loom = loom, key = "Embedding", value = as.data.frame(embedding))
  } else {
    if(sum(c("Embeddings_X","Embeddings_Y")%in%file.h5[["col_attrs"]]$names) != 2) {
      for(i in seq_along(coord.labels)) {
        add_col_attr(loom = loom, key = coord.labels[i], value = as.data.frame(embedding[,i]))
      }
    } else {
      for(i in seq_along(coord.labels)) {
        ca.embeddings<-get_col_attr_by_key(coord.labels[i])
        ca.embeddings<-cbind(ca.embeddings, as.data.frame(embedding[,i]))
        # Update the current coordinates Embeddings
        add_col_attr(loom = loom, key = coord.labels[i], value = as.data.frame(ca.embeddings))
      }
    }
  }
  loom$flush()
  # Ad the given embedding to the global attribute MetaData
  add_global_md_embedding(loom = loom, name = name, is.default = is.default)
}

#########################
# Clusterings functions #
#########################

#'
#'@description Add the given clusters in the given group column attribute and meta data related to the given clustering to the given .loom file handler.
#'@param loom       The loom file handler.
#'@param group      The for the given clustering group to which the given clusters have to be added
#'@param clusters   A list of the the cluster id for each cell present in the matrix
#'
add_clustering<-function(loom
                         , group
                         , clusters) {
  k<-paste0(group,"Clusterings")
  if(col_attrs_exists_by_key(loom = loom, key = k)) {
    ca.clusterings<-get_col_attr_by_key(k)
    clustering<-data.frame(x = clusters)
    colnames(clusters.df)<-paste0("_", ncol(ca.clusterings)+1)
    ca.clusterings<-cbind(ca.clusterings, clustering)
    add_col_attr(loom = loom, key = k, value = as.data.frame(x = ca.clusterings))
  } else {
    clustering<-data.frame("_0" = clusters)
    add_col_attr(loom = loom, key = k, value = as.data.frame(x = clustering))
  }
}

####################
# SCENIC functions #
####################

#'
#'@description Add the regulons with their target genes generated by SCENIC as a row attribute to the given .loom file handler.
#'@param loom     The loom file handler.
#'@param dgem     A matrix of the gene expression with M genes as rows and N cells as columns.
#'@param regulons A list of list of the regulons and their target genes generated by SCENIC.
#'
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
  row.names(reg.mask)<-row.names(dgem)
  colnames(reg.mask)<-gsub(pattern = " ", replacement = "_", x = colnames(reg.mask))
  add_row_attr(loom = loom, key = "Regulons", value = as.data.frame(x = regulons.mask))
}

#'
#'@description Add the regulons AUC matrix generated by SCENIC as a column attribute to the given .loom file handler.
#'@param loom         The loom file handler.
#'@param regulons.AUC A matrix of the regulons AUC values with M regulons as rows and N cells as columns.
#'
add_scenic_regulons_auc_matrix<-function(loom
                                         , regulons.AUC) {
  add_col_attr(loom = loom, key = "RegulonsAUC", value = as.data.frame(x = t(regulons.AUC)))
}

##############################
# Column Meta data functions #
##############################

#'
#'@description Add the Flybase gene as a row attribute to the given .loom file handler.
#'@param loom                       The loom file handler.
#'@param dgem                       A matrix of the gene expression with M genes as rows and N cells as columns.
#'@param fbgn.gn.mapping.file.path  A N-by-2 data.frame containing the mapping between the Flybase gene and the gene symbol.
#'
add_fbgn<-function(loom
                   , dgem
                   , fbgn.gn.mapping.file.path) {
  fbgn.gn.mapping<-read.table(file = fbgn.gn.mapping.file.path, header = F, sep = "\t", quote = '', stringsAsFactors = F)
  colnames(fbgn.gn.mapping)<-c("FBgn","Gene")
  genes<-merge(x = data.frame("Gene"=row.names(dgem)), y = fbgn.gn.mapping, by = "Gene")
  add_row_attr(loom = loom, key = "FBgn", value = genes$FBgn)
}

#####################
# Generic functions #
#####################

#'
#'@description Add a new global attribute to the given .loom object accessible by the given key and containing the given value.
#'@param loom   The loom file handler.
#'@param key    The name of the new added attribute.
#'@param value  The value of the new added attribute.
#'
add_global_attr<-function(loom
                          , key
                          , value) {
  loom$create_attr(attr_name = key, robj = value)
  loom$flush()
}

#'
#'@description Add a new row attribute to the given .loom object accessible by the given key and containing the given value.
#'@param loom   The loom file handler.
#'@param key    The name of the new added attribute.
#'@param value  The value of the new added attribute.
#'
add_row_attr<-function(loom
                       , key
                       , value) {
  ra<-loom[["row_attrs"]]
  ra[[key]]<-value
  loom$flush()
}

col_attrs_exists_by_key<-function(loom
                                , key) {
  ca<-loom[["col_attrs"]]
  return (key %in% names(ca))
}

get_col_attr_by_key<-function(loom
                            , key) {
  ca<-loom[["col_attrs"]]
  return (ca[[key]][])
}

#'
#'@description Add a new column attribute to the given .loom object accessible by the given key and containing the given value.
#'@param loom   The loom file handler.
#'@param key    The name of the new added attribute.
#'@param value  The value of the new added attribute.
#'
add_col_attr<-function(loom
                       , key
                       , value) {
  ca<-loom[["col_attrs"]]
  ca[[key]]<-value
  loom$flush()
}

#'
#'@description Add the given gene expression matrix dgem to the given .loom object.
#'@param loom The loom file handler.
#'@param dgem A matrix of the gene expression with M genes as rows and N cells as columns.
#'
add_matrix<-function(loom
                     , dgem) {
  row.names(dgem)<-NULL
  colnames(dgem)<-NULL
  loom[["matrix"]]<-t(dgem)
  loom$flush()
}

#'
#'@description
#'@loom The loom file handler.
#'
finalize<-function(loom) {
  loom$flush()
  loom$close_all()
}

#'
#'@param file.name                  A string naming the .loom file to be generated.
#'@param dgem                       A matrix of the gene expression with M genes as rows and N cells as columns.
#'@param default.embedding          A M-by-2 data.frame of the embedding (X and Y coordinates) of the cells.
#'@param fbgn.gn.mapping.file.path  A N-by-2 data.frame containing the mapping between the Flybase gene and the gene symbol.
#'
build_loom<-function(file.name
                     , dgem
                     , default.embedding
                     , default.embedding.name
                     , fbgn.gn.mapping.file.path) {
  loom<-H5File$new(filename = file.name, mode = "w")
  cn<-colnames(dgem)
  rn<-row.names(dgem)
  # matrix
  print("Adding matrix...")
  add_matrix(loom = loom, dgem = t(dgem))
  # col_attrs
  print("Adding column attributes...")
  loom$create_group("col_attrs")
  add_col_attr(loom = loom, key = "CellID", value = as.character(cn))
  # Add the default embedding
  add_embedding(loom = loom, embedding = as.data.frame(default.embedding), name = default.embedding.name, is.default = T)
  # row_attrs
  print("Adding row attributes...")
  loom$create_group("row_attrs")
  add_row_attr(loom = loom, key = "Gene", value = as.character(rn))
  add_fbgn(loom = loom, dgem = dgem, fbgn.gn.mapping.file.path = fbgn.gn.mapping.file.path)
  # col_edges
  print("Adding columns edges...")
  col.edges<-loom$create_group("col_edges")
  # row_edges
  print("Adding row edges...")
  row.edges<-loom$create_group("row_edges")
  # layers
  print("Adding layers...")
  layers<-loom$create_group("layers")
  loom$flush()
}
