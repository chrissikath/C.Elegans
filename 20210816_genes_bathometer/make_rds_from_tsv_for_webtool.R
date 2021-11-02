
## Define datasets to work with
#matrices=c(list.files(path=".", pattern=c(".*_data.tsv$",".*_anno.tsv$") ))
matrices=c(list.files(path=".", pattern=c(".*.tsv$") ))

for(matrix in matrices){
    print(matrix)
    data_matrix=as.data.frame(read.table(matrix,header=TRUE,sep="\t",na.strings=c("NA","unknown","","#N/A"),row.names=1))
    print(head(data_matrix))
    ## if you need to transpose the dataframe, the names will be lost
    ## store you data_matrix in df (dataframe)
    #df=data_matrix
    ## Transpose everything, create dfT (T for transposed)
    #dfT <- as.data.frame(as.matrix(t(df)))
    #colnames(dfT)=rownames(data_matrix)
    #data_matrix=dfT
    
    print("Check your matrix' head visually")
    print(data_matrix[1:5,])
    print("Check names(data_matrix) - those should be given, just like the colnames() - these shall not be NULL")
    print(names(data_matrix))

    ## edit the colnames of the data matrix
    #colnames(data_matrix)=gsub("IPI_punch_tumor_","",colnames(data_matrix))
    #colnames(data_matrix)=gsub("_kall_out_b0","",colnames(data_matrix))
    #colnames(data_matrix)=gsub("_re","",colnames(data_matrix))

    ## edit entries in certain columns of the data matrix
    #data_matrix[,"cell.types"]=gsub("\\?","unknown",data_matrix[,"cell.types"])
    #data_matrix$Sex=gsub("Male","M",data_matrix$Sex)
    #data_matrix$Sex=gsub("Female","F",data_matrix$Sex)

    ## edit the rownames of the dtaa matrix
    #new_rn=paste(rownames(data_matrix),data_matrix$gene_name,sep="_")
    #data_matrix[,c(8,9)]=NULL
    
    ## calculate some relation between genes
    #data_matrix["A_over_B",]=data_matrix["A",]/data_matrix["B",]
    
    ## write as .rds
    saveRDS(data_matrix, file=paste0(gsub(".tsv","",matrix),".rds"), ascii = FALSE, version = NULL, compress = TRUE, refhook = NULL)
}

