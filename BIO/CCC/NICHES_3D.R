#改写RunNiches
compute_edgelist <- function (sys.small, position.x, position.y, position.z, k = 6, rad.set = NULL)
{
    df <- data.frame(x = sys.small[[position.x]], y = sys.small[[position.y]], z = sys.small[[position.z]])
    df$barcode <- rownames(df)
    df$x <- as.character(df$x)
    df$y <- as.character(df$y)
    df$z <- as.character(df$z)
    df$x <- as.numeric(df$x)
    df$y <- as.numeric(df$y)
    df$z <- as.numeric(df$z)
    df <- df[, c("x", "y", "z")]
    distance_mat <- apply(df, 1, function(pt) (sqrt(abs(pt["x"] -
        df$x)^2 + abs(pt["y"] - df$y)^2 + abs(pt["z"] - df$z)^2)))
    rownames(distance_mat) <- colnames(distance_mat)
    if (!is.null(k)) {
        message("Compute edgelist based on mutual nearest neighbors.")
        if (!is.null(rad.set))
            warning("'k' is not NULL. Parameter 'rad.set' will be ignored.")
        neighbor_mat <- apply(distance_mat, 1, function(dis_vec) {
            order(dis_vec)[1:(k + 1)]
        })
        adj_mat <- matrix(data = 0, nrow = nrow(distance_mat),
            ncol = ncol(distance_mat))
        rownames(adj_mat) <- rownames(distance_mat)
        colnames(adj_mat) <- colnames(distance_mat)
        for (cell in colnames(neighbor_mat)) adj_mat[cell, neighbor_mat[,
            cell]] <- 1
        adj_mat_final <- 1 * (adj_mat & t(adj_mat)) #mutual nearest neighbors
    }
    else if (!is.null(rad.set)) {
        message("\n Compute edgelist based on spatial radius threshold.")
        adj_mat_final <- 1 * (distance_mat <= rad.set)
    }
    else {
        stop("Both k and rad.set are NULL.")
    }
    edgelist <- igraph::graph.adjacency(adj_mat_final)
    edgelist <- igraph::get.data.frame(edgelist)
    return(edgelist)
}



#' @param object Seurat object. The active identity will be used to define populations for connectomic sampling and crossings.
#' @param assay string. Default: "RNA". The assay to run the NICHES transformation on. 
#' @param LR.database string. Default: "fantom5". Currently accepts "fantom5","omnipath", or "custom".
#' @param species string. The species of the object that is being processed. Only required when LR.database = 'fantom5' with species being 'human','mouse','rat', or 'pig', or LR.database = 'omnipath' with species being 'human','mouse', or 'rat'.
#' @param min.cells.per.ident integer. Default: NULL. A limit on how small (how many cells) a single population can be to participate in connectomic crossings.
#' @param min.cells.per.gene integer. Default: NULL. Limits analysis to interactions involving genes expressed above minimum threshold number of cells in the system. 
#' @param meta.data.to.map character vector. Optional. Default: NULL. A vector of metadata names present in the original object which will be carried to the NICHES objects
#' @param position.x string. Optional. Default: NULL. Only required for spatial omics data. The name that specifies location on the spatial x-axis in the corresponding meta.data column of `object`.
#' @param position.y string. Optional. Default: NULL. Only required for spatial omics data. The name that specifies location on the spatial y-axis in the corresponding meta.data column of `object`.
#' @param custom_LR_database data.frame. Optional. Default: NULL. Only required when LR.database = "custom". Each row is a ligand-receptor mechanism where the first column corresponds to the source genes that express the ligands subunits (separated by '_') and the second column corresponds to the receptor genes that express the receptor subunits (separated by '_').
#' @param k integer. Optional. Default: 4. Number of neighbors in a knn graph. Used to compute a mutual nearest neighbor graph based on the spatial coordinates of the spatial transcriptomic datasets.  
#' @param rad.set numeric. Optional. Default: NULL. The radius threshold to define neighbors based on the spatial coordinates of the spatial transcriptomic datasets. Ignored when 'k' is provided.
#' @param blend string. Default: "mean". Choice of linear operator to combine edges in single-cell niche investigations. Defaults to "mean", also accepts "sum".
#' @param CellToCell logical. Default: TRUE. Whether to analyze cell-cell interactions without considering spatial coordinates.
#' @param CellToSystem logical. Default: FALSE. Whether to analyze summed signaling output to total system coming from each cell. Does not consider Euclidean coordinates.
#' @param SystemToCell logical. Default: FALSE. Whether to analyze summed signaling input from total system landing on each cell (cellular microenvironment/niche). Does not consider Euclidean coordinates. 
#' @param CellToCellSpatial logical. Default: FALSE. Whether to analyze cell-cell interactions between Euclidean neighbors. Only applicable in spatial datasets.
#' @param CellToNeighborhood logical. Default: FALSE. Whether to analyze summed signaling output to Euclidean neighbors. Only applicable in spatial datasets.
#' @param NeighborhoodToCell logical. Default: FALSE. Whether to analyze summed signaling input from Euclidean neighbors (cellular microenvironment/niche). Only applicable in spatial datasets.
#' @param output_format string. Default: "seurat". Choice of the output format. "seurat" will output a list of seurat objects, "raw" will output a list of lists with raw interaction matrix and compiled metadata


RunNICHES.default <- function(object,
                        assay="RNA",
                        LR.database="custom",
                        species,
                        min.cells.per.ident = NULL,
                        min.cells.per.gene = NULL,
                        meta.data.to.map = NULL,
                        position.x = NULL,
                        position.y = NULL,
                        position.z = NULL,
                        custom_LR_database = NULL,
                        k = 6,
                        rad.set = NULL,
                        blend = 'mean',
                        CellToCell = T,
                        CellToSystem = F,
                        SystemToCell = F,
                        CellToCellSpatial = F,
                        CellToNeighborhood = F,
                        NeighborhoodToCell = F,
                        output_format = "seurat",
                        ...
                        ){
  # TODO: check the parameter validity here, then register the parameters
  # 1: check the data format of the required parameters and the optional parameters (if(!is.null())): int, character, etc.
  # 2. check some of the dependencies of the parameters
  #   (1) If LR.database is 'custom', then 'custom_LR_database' can't be NULL
  #   (2) If CellToCellSpatial,CellToNeighborhood, or NeighborhoodToCell is T, then rad.set, 'position.x' and 'position.y' can't be null
  # 3. check whether other functions have done similar checks
  # 4. think whether to put stop or warning
  
  
  # check lr database, species, custom_LR_database
  if(LR.database %in% c("fantom5","omnipath","custom")){
    if(LR.database == "fantom5"){
      if(!(species %in% c('human','mouse','rat','pig'))) 
        stop(paste0("Unsupported species ",species, " for fantom 5 database, only 'human','mouse','rat',or 'pig' supported."))
    }
    if(LR.database == "omnipath"){
      if(!(species %in% c('human','mouse','rat')))
        stop(paste0("Unsupported species ",species, " for omnipath database, only 'human','mouse',or 'rat' supported."))
    }
    if(LR.database == "custom"){
      # check the format of custom_LR_database
      if(is.null(custom_LR_database)) stop("custom_LR_database is NULL")
      # TODO: check gene names
      message("Custom Ligand Receptor database enabled...")
      message("Checking the format of the custom database...")
      if(!is.data.frame(custom_LR_database)){
        warning("Custom database provided is not in dataframe format.")
        warning("Converting to dataframe format...")
        custom_LR_database <- as.data.frame(custom_LR_database)
      }
      if(ncol(custom_LR_database) < 2) stop("Custom database provided contains less than 2 columns.")
    }
  }else stop('\n LR.receptor argument not recognized. Only accepts "omnipath","fantom5" or "custom".')
  
  if(!is.null(custom_LR_database) & LR.database!="custom"){ 
    warning("custom_LR_database is provided but LR.databse is not specified as 'custom'")
}
  # Convert any non-integer inputs to integers. Still allows NULL as option.
  if (!is.null(min.cells.per.ident)){
  min.cells.per.ident <- as.integer(min.cells.per.ident)
  }
  if (!is.null(min.cells.per.gene)){
  min.cells.per.gene <- as.integer(min.cells.per.gene)
  }
  
  # check indicators
  # jc: Add organization names to the list
  org_names_indicator <- c(CellToCell,CellToSystem,SystemToCell,CellToCellSpatial,CellToNeighborhood,NeighborhoodToCell)
  org_names_indicator <- sapply(org_names_indicator,function(org){
    org_out <- as.logical(org)
    if(is.na(org_out)) warning("Organization indicator ",org," is not set to TRUE or FALSE.")
    return(org_out)
  })
  names(org_names_indicator) <- c("CellToCell","CellToSystem","SystemToCell","CellToCellSpatial","CellToNeighborhood","NeighborhoodToCell")
  
  
  # If requested, additionally calculate spatially-limited NICHES organizations
  if (org_names_indicator["CellToCellSpatial"] == T | org_names_indicator["CellToNeighborhood"] == T | org_names_indicator["NeighborhoodToCell"] == T){
    
    if (is.null(position.x) | is.null(position.y) | is.null(position.z)){stop("\n Position information not provided. Please specify metadata columns containing x- and y-axis spatial coordinates.")}
    if(!is.null(k)) k <- as.integer(k)
    if(!is.null(rad.set)) rad.set <- as.numeric(rad.set)
  }
  
  if((!is.null(position.x) | !is.null(position.y) | !is.null(position.z)) & org_names_indicator["CellToCellSpatial"] == F & org_names_indicator["CellToNeighborhood"] == F & org_names_indicator["NeighborhoodToCell"] == F)
    warning("Spatial positions are provided but the spatial organization functions: 'CellToCellSpatial','CellToNeighborhood', and 'NeighborhoodToCell' are set to FALSE.")
    
  
  if(org_names_indicator["CellToSystem"] == T | org_names_indicator["SystemToCell"] == T)
    if(!blend %in% c("sum","mean")) stop("blend paramter is not recognized: need to be 'sum' or 'mean" )
    if(blend == "sum") warning("Operator `sum` will be deprecated in the later release.")
  
  # Initialize output structure
  output <- list()
  
  
  # jc: move the shared preprocessing steps here to avoid redundancy and reduce the number of parameters to be passed to other functions

  # NOTE: relies on Idents(object) to be cell types to subset
  sys.small <- prepSeurat(object,assay,min.cells.per.ident,min.cells.per.gene) 
  ground.truth <- lr_load(LR.database,custom_LR_database,species,rownames(sys.small@assays[[assay]]))
  if (org_names_indicator["CellToCellSpatial"] == T | org_names_indicator["CellToNeighborhood"] == T | org_names_indicator["NeighborhoodToCell"] == T){
    ## 1. Move the neighbor graph construction here
    ## 2. Enable a k-nearest-neighbor parameter as an alternative
    edgelist <- compute_edgelist(sys.small,position.x,position.y,position.z,k,rad.set)
  }
  
  # check the output format
  if(!output_format %in% c("seurat","raw"))
    stop(paste0("Unsupported output format: ",output_format,", Currently only 'seurat' and 'raw' are supported."))

  
  # Calculate NICHES organizations without spatial restrictions
  # jc: only pass the processed data to each function
  # NOTE: RunCellToCell relies on Idents(object) to be cell types to subset
  #       Also each RunXXX function needs Idents(object) to build VectorType meta data 

  if (CellToCell == T){output[[length(output)+1]] <- RunCellToCell(sys.small=sys.small,
                                                                   ground.truth=ground.truth,
                                                                   assay = assay,
                                                                   meta.data.to.map = meta.data.to.map,
                                                                   output_format = output_format
                                                                   )}
  if (CellToSystem == T){output[[length(output)+1]] <- RunCellToSystem(sys.small=sys.small,
                                                                       ground.truth=ground.truth,
                                                                       assay = assay,
                                                                       meta.data.to.map = meta.data.to.map,
                                                                       blend = blend,
                                                                       output_format = output_format
                                                                       )}
  if (SystemToCell == T){output[[length(output)+1]] <- RunSystemToCell(sys.small=sys.small,
                                                                       ground.truth=ground.truth,
                                                                       assay = assay,
                                                                       meta.data.to.map = meta.data.to.map,
                                                                       blend = blend,
                                                                       output_format = output_format
                                                                       )}
  
  
  if (CellToCellSpatial == T){output[[length(output)+1]] <- RunCellToCellSpatial(sys.small=sys.small,
                                                                                 ground.truth=ground.truth,
                                                                                 assay = assay,
                                                                                 meta.data.to.map = meta.data.to.map,
                                                                                 edgelist = edgelist,
                                                                                 output_format = output_format
                                                                                 )} #Spatially-limited Cell-Cell vectors
  if (CellToNeighborhood == T){output[[length(output)+1]] <- RunCellToNeighborhood(sys.small=sys.small,
                                                                                   ground.truth=ground.truth,
                                                                                   assay = assay,
                                                                                   meta.data.to.map = meta.data.to.map,
                                                                                   blend = blend,
                                                                                   edgelist = edgelist,
                                                                                   output_format = output_format
                                                                                   )} #Spatially-limited Cell-Neighborhood vectors
  if (NeighborhoodToCell == T){output[[length(output)+1]] <- RunNeighborhoodToCell(sys.small=sys.small,
                                                                                   ground.truth=ground.truth,
                                                                                   assay = assay,
                                                                                   meta.data.to.map = meta.data.to.map,
                                                                                   blend = blend,
                                                                                   edgelist = edgelist,
                                                                                   output_format = output_format
                                                                                   )} #Spatially-limited Neighborhood-Cell vectors (niches)

  # jc: Add organization names to the list
  names(output) <- names(org_names_indicator)[org_names_indicator]
  
  # Compile objects for output
  return(output)
}



RunNICHES.Seurat <- function(object,
                        assay="RNA",
                        LR.database="cutsom",
                        species,
                        min.cells.per.ident = NULL,
                        min.cells.per.gene = NULL,
                        meta.data.to.map = NULL,
                        position.x = NULL,
                        position.y = NULL,
                        position.z = NULL,
                        cell_types = NULL,
                        custom_LR_database = NULL,
                        k = 6,
                        rad.set = NULL,
                        blend = 'mean',
                        CellToCell = T,
                        CellToSystem = F,
                        SystemToCell = F,
                        CellToCellSpatial = F,
                        CellToNeighborhood = F,
                        NeighborhoodToCell = F,
                        output_format = "seurat",
                        ...
                        ){
  # TODO: check the parameter validity here, thenf register the parameters
  # 1: check the data format of the required parameters and the optional parameters (if(!is.null())): int, character, etc.
  # 2. check some of the dependencies of the parameters
  #   (1) If LR.database is 'custom', then 'custom_LR_database' can't be NULL
  #   (2) If CellToCellSpatial,CellToNeighborhood, or NeighborhoodToCell is T, then rad.set, 'position.x' and 'position.y' can't be null
  # 3. check whether other functions have done similar checks
  # 4. think whether to put stop or warning
  
  # check assay:
  if(is.null(object@assays[[assay]])) stop(paste0("Assay ",assay," is NULL"))
  
  
  # check meta.data.to.map
  if(!is.null(meta.data.to.map)){
    # check each meta data is present in the provided seurat object
    for(each_meta in meta.data.to.map){
      if(!each_meta %in% colnames(object@meta.data)) 
        warning(paste0("Metadata: "),each_meta," is not present in the provided Seurat object")
    }
  }

  # Alternative: stop when cell_types not provided
  # check cell types 
  if(is.null(cell_types)){
    if(CellToCell == T) stop("cell_types need to be provided to run the downsampling procedure in the CellToCell computation") 
    else warning("cell_types not provided. The Identity of the object is assumed to be cell types")
  }
  else if(!(cell_types %in% colnames(object@meta.data))){
    if(CellToCell == T) stop("cell_types is not in the meta.data columns of the input Seurat object,please make sure its name is in the colnames of meta.data")
    else warning("cell_types is not in the meta.data columns of the input Seurat object. The Identity of the object is assumed to be cell types")
  }
  else{
    message("Set cell types as Identity of object internally")
    Seurat::Idents(object) <- cell_types
  }
  # check cell types metadata when CellToCell is True
  #if(CellToCell == T){
  #  if(is.null(cell_types)) stop("cell_types need to be provided to run the downsampling procedure in the CellToCell computation") 
  #  else if(!(cell_types %in% colnames(object@meta.data))) stop("cell_types is not in the meta.data columns of the input Seurat object, 
  #                                                           please make sure its name is in the colnames of meta.data")
  #}
  # set the ident to cell_types
  output <- RunNICHES.default(object = object,
                        assay = assay,
                        LR.database = LR.database,
                        species = species,
                        min.cells.per.ident = min.cells.per.ident,
                        min.cells.per.gene = min.cells.per.gene,
                        meta.data.to.map = meta.data.to.map,
                        position.x = position.x,
                        position.y = position.y,
                        position.z = position.z,
                        custom_LR_database = custom_LR_database,
                        k = k,
                        rad.set = rad.set,
                        blend = blend,
                        CellToCell = CellToCell,
                        CellToSystem = CellToSystem,
                        SystemToCell = SystemToCell,
                        CellToCellSpatial = CellToCellSpatial,
                        CellToNeighborhood = CellToNeighborhood,
                        NeighborhoodToCell = NeighborhoodToCell,
                        output_format = output_format,
                        ...)
  # Compile objects for output
  return(output)
}


RunNeighborhoodToCell <- function(sys.small, ground.truth, assay, meta.data.to.map, blend = "mean", edgelist, output_format)
{
    subunit.list <- list()
    for (s in 1:ncol(ground.truth$source.subunits)) {
        subunit.list[[s]] <- matrix(data = 1, nrow = nrow(ground.truth$source.subunits),
            ncol = ncol(sys.small@assays[[assay]]@data[, edgelist$from]))
        colnames(subunit.list[[s]]) <- colnames(sys.small@assays[[assay]]@data[,
            edgelist$from])
        rownames(subunit.list[[s]]) <- rownames(ground.truth$source.subunits)
        non.na.indices <- !is.na(ground.truth$source.subunits[,
            s])
        subunit.list[[s]][non.na.indices, ] <- as.matrix(sys.small@assays[[assay]]@data[ground.truth$source.subunits[non.na.indices,
            s], edgelist$from])
    }
    lig.data <- Reduce("*", subunit.list)
    rm(subunit.list)
    subunit.list <- list()
    for (t in 1:ncol(ground.truth$target.subunits)) {
        subunit.list[[t]] <- matrix(data = 1, nrow = nrow(ground.truth$target.subunits),
            ncol = ncol(sys.small@assays[[assay]]@data[, edgelist$to]))
        colnames(subunit.list[[t]]) <- colnames(sys.small@assays[[assay]]@data[,
            edgelist$to])
        rownames(subunit.list[[t]]) <- rownames(ground.truth$target.subunits)
        non.na.indices <- !is.na(ground.truth$target.subunits[,
            t])
        subunit.list[[t]][non.na.indices, ] <- as.matrix(sys.small@assays[[assay]]@data[ground.truth$target.subunits[non.na.indices,
            t], edgelist$to])
    }
    rec.data <- Reduce("*", subunit.list)
    rm(subunit.list)
    scc <- lig.data * rec.data
    rownames(scc) <- paste(rownames(lig.data), rownames(rec.data),
        sep = "—")
    colnames(scc) <- colnames(rec.data)
    scc <- as.matrix(scc)
    if (blend == "sum")
        scc <- t(rowsum(t(scc), colnames(scc)))
    else if (blend == "mean")
        scc <- sapply(unique(colnames(scc)), function(receiving_name) rowMeans(scc[,
            colnames(scc) == receiving_name, drop = FALSE], na.rm = TRUE))
    barcodes <- colnames(scc)
    colnames(scc) <- paste("Neighborhood", colnames(scc), sep = "—")
    demo <- Seurat::CreateSeuratObject(counts = as.matrix(scc),
        assay = "NeighborhoodToCell")
    demo <- Seurat::AddMetaData(demo, metadata = barcodes, col.name = "ReceivingCell")
    demo <- Seurat::AddMetaData(demo, metadata = Seurat::Idents(sys.small)[barcodes],
        col.name = "ReceivingType")
    if (!is.null(meta.data.to.map)) {
        receiving.barcodes <- barcodes
        receiving.metadata <- as.matrix(sys.small@meta.data[,
            meta.data.to.map, drop = FALSE][receiving.barcodes,
            ])
        meta.data.to.add.also <- receiving.metadata
        rownames(meta.data.to.add.also) <- paste("Neighborhood",
            receiving.barcodes, sep = "—")
        demo <- Seurat::AddMetaData(demo, metadata = as.data.frame(meta.data.to.add.also))
    }
    Seurat::Idents(demo) <- demo$ReceivingType
    message(paste("\n", length(unique(demo$ReceivingCell)), "Neighborhood-To-Cell edges were computed, across",
        length(unique(demo$ReceivingType)), "cell types"))
    if (output_format == "seurat")
        return(demo)
    else {
        output_list <- vector(mode = "list", length = 2)
        names(output_list) <- c("NeighborhoodToCellMatrix", "metadata")
        output_list[["NeighborhoodToCellMatrix"]] <- demo[["NeighborhoodToCell"]]@counts
        output_list[["metadata"]] <- demo@meta.data
        return(output_list)
    }
}



RunCellToCellSpatial <- function(sys.small,
                                 ground.truth,
                                 assay,
                                 meta.data.to.map,
                                 edgelist,
                                 output_format
                                 ){

  
    # Make ligand matrix

    #lig.data <- sys.small@assays[[assay]]@data[ligands,edgelist$from]

    subunit.list <- list() # Builds sending (ligand) data for any number of ligand subunits
    for (s in 1:ncol(ground.truth$source.subunits)){ #For each subunit column...
      subunit.list[[s]] <- matrix(data = 1,nrow = nrow(ground.truth$source.subunits),ncol = ncol(sys.small@assays[[assay]]@data[,edgelist$from])) #initialize a mechanism x barcode matrix of all NAs
      colnames(subunit.list[[s]]) <- colnames(sys.small@assays[[assay]]@data[,edgelist$from])
      rownames(subunit.list[[s]]) <- rownames(ground.truth$source.subunits)
      non.na.indices <- !is.na(ground.truth$source.subunits[,s]) #Identify rows in the s-th column of the ground truth which are not NA
      subunit.list[[s]][non.na.indices,] <- as.matrix(sys.small@assays[[assay]]@data[ground.truth$source.subunits[non.na.indices,s],edgelist$from])   #For every row in the initialized matrix corresponding to the indices of the ground.truth which are not NA, replace with the rows from the Seurat object corresponding to the genes in the ground.truth at those indices
    }
    lig.data <- Reduce('*',subunit.list)
    rm(subunit.list)
    
  # Make receptor matrix

    #rec.data <- sys.small@assays[[assay]]@data[receptors,edgelist$to]
    
    subunit.list <- list() # Builds receiving (receptor) data for any number of receptor subunits
    for (t in 1:ncol(ground.truth$target.subunits)){
      subunit.list[[t]] <- matrix(data = 1,nrow = nrow(ground.truth$target.subunits),ncol = ncol(sys.small@assays[[assay]]@data[,edgelist$to])) #initialize a mechanism x barcode matrix of all NAs
      colnames(subunit.list[[t]]) <- colnames(sys.small@assays[[assay]]@data[,edgelist$to])
      rownames(subunit.list[[t]]) <- rownames(ground.truth$target.subunits)
      non.na.indices <- !is.na(ground.truth$target.subunits[,t]) #Identify rows in the t-th column of the ground truth which are not NA
      subunit.list[[t]][non.na.indices,] <- as.matrix(sys.small@assays[[assay]]@data[ground.truth$target.subunits[non.na.indices,t],edgelist$to])   #For every row in the initialized matrix corresponding to the indices of the ground.truth which are not NA, replace with the rows from the Seurat object corresponding to the genes in the ground.truth at those indices
    }
    rec.data <- Reduce('*',subunit.list)
    rm(subunit.list)

  # Make SCC matrix
  scc <- lig.data*rec.data
  rownames(scc) <- paste(rownames(lig.data),rownames(rec.data),sep = '—')
  colnames(scc) <- paste(colnames(lig.data),colnames(rec.data),sep = '—')
  sending.cell.idents <- as.character(Seurat::Idents(sys.small)[colnames(lig.data)])
  receiving.cell.idents <- as.character(Seurat::Idents(sys.small)[colnames(rec.data)])
  dim(scc)

  # Use this matrix to create a Seurat object:
  demo <- Seurat::CreateSeuratObject(counts = as.matrix(scc),assay = 'CellToCellSpatial')

  # Add key metadata

  meta.data.to.add <- data.frame(SendingType = sending.cell.idents,
                                 ReceivingType = receiving.cell.idents)

  rownames(meta.data.to.add) <- colnames(scc)
  meta.data.to.add$VectorType <- paste(meta.data.to.add$SendingType,
                                       meta.data.to.add$ReceivingType,
                                       sep = '—')

  #Add metadata to the Seurat object
  demo <- Seurat::AddMetaData(demo,metadata = meta.data.to.add)

  # Gather and assemble additional metadata
  if (!is.null(meta.data.to.map)){
    # Identify sending and receiving barcodes
    sending.barcodes <- colnames(lig.data)
    receiving.barcodes <- colnames(rec.data)
    # Pull and format sending and receiving metadata
    # jc: possible bug, change object to sys.small
    sending.metadata <- as.matrix(sys.small@meta.data[,meta.data.to.map,drop=FALSE][sending.barcodes,])
    receiving.metadata <- as.matrix(sys.small@meta.data[,meta.data.to.map,drop=FALSE][receiving.barcodes,])
    # Make joint metadata
    datArray <- abind::abind(sending.metadata,receiving.metadata,along=3)
    joint.metadata <- as.matrix(apply(datArray,1:2,function(x)paste(x[1],"-",x[2])))
    # Define column names
    colnames(joint.metadata) <- paste(colnames(sending.metadata),'Joint',sep = '.')
    colnames(sending.metadata) <- paste(colnames(sending.metadata),'Sending',sep='.')
    colnames(receiving.metadata) <- paste(colnames(receiving.metadata),'Receiving',sep='.')
    # Compile
    meta.data.to.add.also <- cbind(sending.metadata,receiving.metadata,joint.metadata)
    rownames(meta.data.to.add.also) <- paste(sending.barcodes,receiving.barcodes,sep='—')
    # Add additional metadata
    demo <- Seurat::AddMetaData(demo,metadata = as.data.frame(meta.data.to.add.also))
  }
  # Set initial identity
  Seurat::Idents(demo) <- demo$VectorType
  
  # How many vectors were captured by this sampling?
  message(paste("\n",length(unique(demo$VectorType)),'distinct VectorTypes were computed, out of',length(table(Seurat::Idents(sys.small)))^2,'total possible'))

  
  if(output_format == "seurat") return(demo)
  else{
    output_list <- vector(mode = "list",length=2)
    names(output_list) <- c("CellToCellSpatialMatrix","metadata")
    output_list[["CellToCellSpatialMatrix"]] <- demo[["CellToCellSpatial"]]@counts
    output_list[["metadata"]] <- demo@meta.data
    return(output_list)
  }
  
}

RunCellToNeighborhood <- function(sys.small,
                                  ground.truth,
                                  assay,
                                  meta.data.to.map,
                                  blend="mean",
                                  edgelist,
                                  output_format
                                  ){

  # Make ligand matrix
  
  #lig.data <- sys.small@assays[[assay]]@data[ligands,edgelist$from]
  
  subunit.list <- list() # Builds sending (ligand) data for any number of ligand subunits
  for (s in 1:ncol(ground.truth$source.subunits)){ #For each subunit column...
    subunit.list[[s]] <- matrix(data = 1,nrow = nrow(ground.truth$source.subunits),ncol = ncol(sys.small@assays[[assay]]@data[,edgelist$from])) #initialize a mechanism x barcode matrix of all NAs
    colnames(subunit.list[[s]]) <- colnames(sys.small@assays[[assay]]@data[,edgelist$from])
    rownames(subunit.list[[s]]) <- rownames(ground.truth$source.subunits)
    non.na.indices <- !is.na(ground.truth$source.subunits[,s]) #Identify rows in the s-th column of the ground truth which are not NA
    subunit.list[[s]][non.na.indices,] <- as.matrix(sys.small@assays[[assay]]@data[ground.truth$source.subunits[non.na.indices,s],edgelist$from])   #For every row in the initialized matrix corresponding to the indices of the ground.truth which are not NA, replace with the rows from the Seurat object corresponding to the genes in the ground.truth at those indices
  }
  lig.data <- Reduce('*',subunit.list)
  rm(subunit.list)
  
  # Make receptor matrix
  
  #rec.data <- sys.small@assays[[assay]]@data[receptors,edgelist$to]
  
  subunit.list <- list() # Builds receiving (receptor) data for any number of receptor subunits
  for (t in 1:ncol(ground.truth$target.subunits)){
    subunit.list[[t]] <- matrix(data = 1,nrow = nrow(ground.truth$target.subunits),ncol = ncol(sys.small@assays[[assay]]@data[,edgelist$to])) #initialize a mechanism x barcode matrix of all NAs
    colnames(subunit.list[[t]]) <- colnames(sys.small@assays[[assay]]@data[,edgelist$to])
    rownames(subunit.list[[t]]) <- rownames(ground.truth$target.subunits)
    non.na.indices <- !is.na(ground.truth$target.subunits[,t]) #Identify rows in the t-th column of the ground truth which are not NA
    subunit.list[[t]][non.na.indices,] <- as.matrix(sys.small@assays[[assay]]@data[ground.truth$target.subunits[non.na.indices,t],edgelist$to])   #For every row in the initialized matrix corresponding to the indices of the ground.truth which are not NA, replace with the rows from the Seurat object corresponding to the genes in the ground.truth at those indices
  }
  rec.data <- Reduce('*',subunit.list)
  rm(subunit.list)
  
  
  # Make SCC matrix
  scc <- lig.data*rec.data
  rownames(scc) <- paste(rownames(lig.data),rownames(rec.data),sep = '—')

  # Condense by column name
  colnames(scc) <- colnames(lig.data) # Make colnames equal to sending cell
  scc <- as.matrix(scc)
  if(blend == "sum") scc <- t(rowsum(t(scc), colnames(scc)))
  else if(blend == "mean")  scc <- sapply(unique(colnames(scc)), function(sending_name) 
    rowMeans(scc[,colnames(scc)== sending_name,drop=FALSE], na.rm=TRUE) )
  
  # Label columns properly
  barcodes <- colnames(scc)
  colnames(scc) <- paste(colnames(scc),'Neighborhood',sep = '—')

  # Use this matrix to create a Seurat object:
  demo <- Seurat::CreateSeuratObject(counts = as.matrix(scc),assay = 'CellToNeighborhood')

  # Add metadata based on ident slot
  demo <- Seurat::AddMetaData(demo,metadata = barcodes,col.name = 'SendingCell')
  demo <- Seurat::AddMetaData(demo,metadata = Seurat::Idents(sys.small)[barcodes],col.name = 'SendingType')

  # Gather and assemble additional metadata
  if (!is.null(meta.data.to.map)){
    # Identify sending and receiving barcodes
    sending.barcodes <- barcodes # Only sending cell metadata applies for this function
    #receiving.barcodes <- colnames(rec.map)
    # Pull and format sending and receiving metadata
    # jc: possible bug, change object to sys.small
    sending.metadata <- as.matrix(sys.small@meta.data[,meta.data.to.map,drop=FALSE][sending.barcodes,])
    #receiving.metadata <- as.matrix(object@meta.data[,meta.data.to.map][receiving.barcodes,])
    # Make joint metadata
    #datArray <- abind(sending.metadata,receiving.metadata,along=3)
    #joint.metadata <- as.matrix(apply(datArray,1:2,function(x)paste(x[1],"-",x[2])))
    # Define column names
    #colnames(joint.metadata) <- paste(colnames(sending.metadata),'Joint',sep = '.')
    #colnames(sending.metadata) <- paste(colnames(sending.metadata),'Sending',sep='.')
    #colnames(receiving.metadata) <- paste(colnames(receiving.metadata),'Receiving',sep='.')
    # Compile
    meta.data.to.add.also <- sending.metadata
    rownames(meta.data.to.add.also) <- paste(sending.barcodes,'Neighborhood',sep='—')
    # Add additional metadata
    demo <- Seurat::AddMetaData(demo,metadata = as.data.frame(meta.data.to.add.also))
  }
  # Set initial identity
  Seurat::Idents(demo) <- demo$SendingType
  # How many vectors were captured by this sampling?
  message(paste("\n",length(unique(demo$SendingCell)),'Cell-To-Neighborhood edges were computed, across',length(unique(demo$SendingType)),'cell types'))

  if(output_format == "seurat") return(demo)
  else{
    output_list <- vector(mode = "list",length=2)
    names(output_list) <- c("CellToNeighborhoodMatrix","metadata")
    output_list[["CellToNeighborhoodMatrix"]] <- demo[["CellToNeighborhood"]]@counts
    output_list[["metadata"]] <- demo@meta.data
    return(output_list)
  }
  
}