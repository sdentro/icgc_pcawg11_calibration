morris_file = "data/morris/purity_ploidy.txt"
vanloo_wedge_file = "data/vanloo_wedge/1_purity_ploidy/purity_ploidy.txt"

morris = parse.purity.ploidy(morris_file)
vanloo_wedge = parse.purity.ploidy(vanloo_wedge_file)
list_of_tables = list(morris, vanloo_wedge)
vector_of_names = c("morris", "vanloo_wedge")

#' Convert a list of tables into a matrix with purity estimates only, placing NAs where a group has not reported anything
create.purity.table = function(list_of_tables, vector_of_names) {
  purity = list_of_tables[[1]][,c("sample", "purity")]
  for (i in 2:length(list_of_tables)) {
    # Add a new column for this purity table
    new.col = list_of_tables[[i]]$purity[match(purity$sample, list_of_tables[[i]]$sample)]
    purity = cbind(purity, new.col)
    
    # Create and add rows for samples not previously mentioned
    new.index = which(!(list_of_tables[[i]]$sample %in% purity$sample))
    new.rows = matrix(NA, ncol=ncol(purity), nrow=length(new.index))
    new.rows[,1] = list_of_tables[[i]]$sample[new.index]
    new.rows[,i+1] = list_of_tables[[i]]$purity[new.index]
    colnames(new.rows) = colnames(purity)
    purity = rbind(purity, new.rows)
  }
  colnames(purity) = c("sample", vector_of_names)
  return(purity[order(purity$sample),])
}

#' Create a ploidy overview table, inserting NAs where ploidy was not reported
create.ploidy.table = function(list_of_tables, vector_of_names) {
  ploidy = list_of_tables[[1]][,c("sample", "ploidy")]
  for (i in 2:length(list_of_tables)) {
    # Add a new column for this ploidy table
    new.col = list_of_tables[[i]]$ploidy[match(ploidy$sample, list_of_tables[[i]]$sample)]
    ploidy = cbind(ploidy, new.col)
    
    # Create and add rows for samples not previously mentioned
    new.index = which(!(list_of_tables[[i]]$sample %in% ploidy$sample))
    new.rows = matrix(NA, ncol=ncol(ploidy), nrow=length(new.index))
    new.rows[,1] = list_of_tables[[i]]$sample[new.index]
    new.rows[,i+1] = list_of_tables[[i]]$ploidy[new.index]
    colnames(new.rows) = colnames(ploidy)
    ploidy = rbind(ploidy, new.rows)
  }
  colnames(ploidy) = c("sample", vector_of_names)
  return(ploidy[order(ploidy$sample),])
}

purity = create.purity.table(list_of_tables, vector_of_names)
write.table(purity, "1_purity_ploidy/purity.tsv", sep="\t", quote=F, row.names=F)
ploidy = create.ploidy.table(list_of_tables, vector_of_names)
write.table(ploidy, "1_purity_ploidy/ploidy.tsv", sep="\t", quote=F, row.names=F)