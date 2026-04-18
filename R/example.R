load("/work/tingthou/SingleCell/RealData/artista_data.RData")
make_celltype_mat = function(celltype_annotation){

  celltype_mat = matrix(0,length(celltype_annotation),length(c(table(celltype_annotation))))

  for(i in 1:length(table(celltype_annotation))){
    celltypename = names(table(celltype_annotation))[i]
    celltype_mat[,i] = 0
    celltype_mat[which(celltype_annotation==celltypename),i] =1
  }
  colnames(celltype_mat) = names(table(celltype_annotation))
  return(celltype_mat)

}





# Now we can create the cell type matrix, just like the cell type proportion matrix we had in the spot level data, but now for single cell level data, the values are binary indicators rather than 0-1 proportion values.
data_use$celltype = make_celltype_mat(data_use$cluster)
rownames(data_use$celltype) = colnames(data_use$count_use)

Obj = Create_SIGNAL_Object(cell_type_compositions = data_use$celltype[1:1000,],
                           gene_expression = as.matrix(data_use$count_use[1:40,1:1000]),
                           location = as.matrix(data_use$location_coord[1:1000,]),
                           covariates = NULL,
                           project = "Artista")

Obj<-data_preprocess(Obj)
Obj = Calculate_Kernel(Obj, approximation = FALSE)
Obj = Calculate_celltype_Kernel(Obj)
a<-RunTesting(Obj,num_cores = 30)
