source("Parser.R")
source("Plotting.R")

samplename = "0c7af04b-e171-47c4-8be5-5db33f20148e"
morris_file = "data/morris/mutation_assignment/mutation_assignment.all.0c7af04b-e171-47c4-8be5-5db33f20148e.csv"
vanloo_wedge_file = "data/vanloo_wedge/2_clustering/0c7af04b-e171-47c4-8be5-5db33f20148e/0c7af04b-e171-47c4-8be5-5db33f20148e_cluster_membership.txt"

morris = parse.mut.assignments(morris_file)
vanloo_wedge = parse.mut.assignments(vanloo_wedge_file)

vector_of_names = c("morris", "vanloo_wedge")
list_of_tables = list(morris, vanloo_wedge)

#######################################################################
# Get which mutations were assigned by whome
#######################################################################
get.mut.assigned.by.who = function(master_table, list_of_tables, vector_of_names) {
  for (i in 1:length(list_of_tables)) {
    #is_assigned = array(F, nrow(master_table))
    sub = list_of_tables[[i]]
    num.cols = ncol(sub)
    # Subset such that all mentioned mutations are assigned
    is_assigned = sapply(1:nrow(sub), function(i, sub, num.cols) { !any(is.na(sub[i,3:num.cols])) }, sub=sub, num.cols=num.cols)
    sub = sub[is_assigned,]
    master_table = cbind(master_table, (master_table[,1] %in% sub[,1] & master_table[,2] %in% sub[,2]))
  }
  colnames(master_table) = c("Chromosome", "Position", vector_of_names)
  return(master_table)
}

master_table = vanloo_wedge[,c(1,2)]
assignment_overview = get.mut.assigned.by.who(master_table, list_of_tables, vector_of_names)
write.table(assignment_overview, paste("2_mutation_assignments/tables/", samplename, "_mutation_assignment_inventory.tsv", sep=""), sep="\t", quote=F, row.names=F)

#######################################################################
# Synchronise mutation assignment tables
#######################################################################
#' Tests which mutations were assigned by all methods.
#' @return A vector with a boolean for each mutation.
mutation.assigned.by.all = function(assignment.tables, raw.data, vector_of_names) {
  sel = rep(T, nrow(raw.data))
  for (i in 1:length(vector_of_names)) {
    # Don't take into account empty assignment tables, which occurs when a method hasn't pruduced output for a sample
    if (nrow(assignment.tables[[i]]) > 0) {
      sel = sel & (raw.data$Chromosome %in% assignment.tables[[i]][,1] & raw.data$Position %in% assignment.tables[[i]][,2])
    }
  }
  return(sel)
}

#' Synchronise the assignments based on whether a mutation is assigned by all methods. Also annotates Subclonal.fraction to the df.
sync.assignments = function(assignment.tables, raw.data, selection, vector_of_names) {
  dat.shared = list()
  for (i in 1:length(vector_of_names)) {
    dat = assignment.tables[[i]]
    if (nrow(dat) == 0) {
      # This happens when a method hasn't produced results for a sample. Reuse the empty data.frame
      dat.shared[[i]] = dat
    } else {
      dat.sel = dat[,1] %in% raw.data$Chromosome[selection] & dat[,2] %in% raw.data$Position[selection]
      dat.shared[[i]] = dat[dat.sel,]
      dat.shared[[i]]$Subclonal.fraction = raw.data$Subclonal.fraction[selection]
    }
  }
  return(dat.shared)
}

#' Calculate the identity matrix for all methods
calc.ident.matrices = function(dat.shared, vector_of_names) {
  ident.matrices = list()
  for (i in 1:length(vector_of_names)) {
    if (nrow(dat.shared[[i]]) == 0) {
      # If there are no mutation assignments available, set the identity matrix to NA
      ident.matrices[[i]] = NULL
    } else {
      #ident.matrices[[i]] = GetIdentityArrayFromAssignments(dat.shared[[i]]$Cluster)
      ident.matrices[[i]] = GetIdentityArrayFromProbabilities(as.matrix(dat.shared[[i]][,3:(ncol(dat.shared[[i]])-1)]))
    }
  }
  return(ident.matrices)
}

d = parse.dp.input(paste("data/dirichlet_input/", samplename, "_allDirichletProcessInfo.txt", sep=""))

# Sync the mutation calls between tables
# Only take mutations that are used by all methods
sel = mutation.assigned.by.all(list_of_tables, d, vector_of_names)

# Subset the assignment tables using the above selection to get the shared mutations
dat.shared = sync.assignments(list_of_tables, d, sel, vector_of_names)

# Calculate the identity matrix
ident.matrices = calc.ident.matrices(dat.shared, vector_of_names)

#######################################################################
# Calculate similarities between each pair of matrices
#######################################################################
diff.mse = matrix(0, length(vector_of_names), length(vector_of_names))
for (i in 1:length(vector_of_names)) {
  for (j in i:length(vector_of_names)) {
    if (!is.null(ident.matrices[[i]]) & !is.null(ident.matrices[[j]])) {
      diff.mse[i, j] = mean((ident.matrices[[i]] - ident.matrices[[j]])^2)
    } else {
      diff.mse[i, j] = NA
    }
  }
}
diff.mse.d = data.frame(diff.mse)
colnames(diff.mse.d) = vector_of_names
write.table(diff.mse.d, file=paste("2_mutation_assignments/similarities/", samplename, ".mse.txt", sep=""), sep="\t", row.names=F, quote=F)


#######################################################################
# Plot a heatmap with data ordered according to the first method
#######################################################################
most.likely.node.assignments = apply(dat.shared[[1]][,3:(ncol(dat.shared[[1]])-1)], 1, which.max)
ord = order(most.likely.node.assignments)
for (i in 1:length(vector_of_names)) {
  # Only plot heatmap when there is a matrix, i.e. when a method has produced mutation assignments for this sample
  if (!is.null(ident.matrices[[i]])) {
    if (nrow(ident.matrices[[i]]) < 10000) {
      plotHeatmapFull(ident.matrices[[i]][ord,ord], paste("figures/", samplename, "_", vector_of_names[i], ".png", sep=""))
    } else {
      plotHeatmapMedium(ident.matrices[[i]][ord,ord], paste("figures/", samplename, "_", vector_of_names[i], ".png", sep=""))
    }
  }
}

#######################################################################
# Save a table with all most likely assignments
#######################################################################
most.likely.node.assignments = cbind(dat.shared[[1]][,c(1,2)], most.likely.node.assignments)
for (i in 2:length(dat.shared)) {
  most.likely.node.assignments = cbind(most.likely.node.assignments, apply(dat.shared[[i]][,3:(ncol(dat.shared[[i]])-1)], 1, which.max))
}
colnames(most.likely.node.assignments) = c("chr", "pos", vector_of_names)
write.table(most.likely.node.assignments, file=paste("2_mutation_assignments/tables/", samplename, "_most_likely_node_assignments.tsv", sep=""), sep="\t", col.names=T, quote=F, row.names=F)

#######################################################################
# Obtain cluster locations by taking the mean Subclonal fraction of its assigned mutations
#######################################################################
#' Get mean cluster locations based on the Subclonal fraction estimates of its assigned mutations
get.cluster.locations.shared.muts = function(dat.shared, vector_of_names) {
  mean.cluster.locs = list()
  for (i in 1:length(vector_of_names)) {
    dd = dat.shared[[i]]
    if (nrow(dd) == 0) {
      # No clusters available when a method hasn't produced mutation assignments.
      mean.cluster.locs[[i]] = NULL
    } else {
      assignment = apply(dd[,3:(ncol(dd)-1)], 1, which.max)
      mean.subcl.frac = data.frame()
      for (cluster in unique(assignment)) {
        cluster.subcl.frac = dd$Subclonal.fraction[assignment==cluster]
        mean.subcl.frac = rbind(mean.subcl.frac, c(cluster, mean(cluster.subcl.frac), mean(cluster.subcl.frac)+sd(cluster.subcl.frac), sum(assignment==cluster)))
      }
      colnames(mean.subcl.frac) = c("Cluster", "Subclonal.fraction.mean", "Subclonal.fraction.sd", "Total.muts")
      mean.cluster.locs[[i]] = mean.subcl.frac
    }
  }
  return(mean.cluster.locs)
}

cluster.locations = get.cluster.locations.shared.muts(dat.shared, vector_of_names)
for (i in 1:length(vector_of_names)) {
  if (!is.null(cluster.locations[[i]])) {
    write.table(file=paste("2_mutation_assignments/tables/", samplename, "_cluster_locations_shared_", vector_of_names[i], ".txt", sep=""), cluster.locations[[i]], row.names=F, quote=F, sep="\t")
  }
}

#' Obtain average CCF of clusters using all assigned mutations
get.cluster.locations.all.muts = function(list_of_tables, raw.data, vector_of_names) {
  mean.cluster.locs = list()
  for (i in 1:length(list_of_tables)) {
    dd = list_of_tables[[1]]
    selection = d$Chromosome %in% dd[,1] & d$Position %in% dd[,2]
    dd = cbind(dd, d$Subclonal.fraction[selection])
    assignment = apply(d[,3:(ncol(d)-1)], 1, which.max)
    
    mean.subcl.frac = data.frame()
    for (cluster in unique(assignment)) {
      cluster.subcl.frac = dd[assignment==cluster, ncol(dd)]
      mean.subcl.frac = rbind(mean.subcl.frac, c(cluster, mean(cluster.subcl.frac), mean(cluster.subcl.frac)+sd(cluster.subcl.frac), sum(assignment==cluster)))
    }
    colnames(mean.subcl.frac) = c("Cluster", "Subclonal.fraction.mean", "Subclonal.fraction.sd", "Total.muts")
    mean.cluster.locs[[i]] = mean.subcl.frac
  }
  return(mean.cluster.locs)
}

cluster.locations = get.cluster.locations.all.muts(list_of_tables, d, vector_of_names)
for (i in 1:length(vector_of_names)) {
  if (!is.null(cluster.locations[[i]])) {
    write.table(file=paste("2_mutation_assignments/tables/", samplename, "_cluster_locations_all_", vector_of_names[i], ".txt", sep=""), cluster.locations[[i]], row.names=F, quote=F, sep="\t")
  }
}