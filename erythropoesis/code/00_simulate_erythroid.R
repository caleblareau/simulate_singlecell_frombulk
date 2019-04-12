library(data.table)
library(dplyr)
library(Matrix)
library(BuenColors)
library(stringr)
library(cowplot)
library(irlba)

eryth_color_maps <- c("P1" = "#3b82ae", "P2" = "#547294", "P3" = "#6d617a", "P4" = "#865160", "P5" = "#9f4046", "P6" = "#b8302c", "P7" = "#d11f12", "P8" = "#de1705")
all_color_maps <- c(ejc_color_maps[c("HSC", "MPP", "CMP")], "MEP" = "#FF81AF" ,eryth_color_maps)

#-----------
# Parameters
#-----------
# n_cells Number of cells from each group to be simulated; either a number (all the same) or a vector
# of length(which_celltypes)

# which_celltypes Number of groups to partition cells into (even groups); must evenly divide into n_cells
# n_frags_per_cell number of ragments in peaks to be simulated per single cell

# rate_noise number between 0 (perfect downsample) and 1 (nonsense) for noise

# seed rando parameter for setting the seed
# shuffle Randomly order the resulting cells and peaks

bulk <- data.matrix(data.frame(fread("../data/combined_12.counts.tsv")))

# Reorder based on trajectories
bulk <- bulk[,c(2,4,1,3, 5:12)]

simulate_scatac <- function(n_cells, which_celltypes, n_frags_per_cell = 1000, 
                            rate_noise = 0, seed = 100, shuffle = FALSE){
  
  # Reproducibility
  set.seed(seed)
  which_celltypes <- sort(which_celltypes)
  stopifnot(rate_noise < 1) 
  stopifnot(n_frags_per_cell > 100)
  n_peaks <- dim(bulk)[1]
  #--
  # Set up cell labels
  #--
  
  if(length(n_cells) > 1){
    stopifnot(length(which_celltypes) == length(n_cells))
    
    # Generate cell labels
    cell_labels <- sapply(1:length(which_celltypes), function(i){
      rep(which_celltypes[i], n_cells[i])
    }) %>% unlist() %>% sort()
    
  } else {
    n_groups <- length(which_celltypes)
    cell_labels <- sort(rep(which_celltypes, n_cells*n_groups))
  }
  final_names <- paste0(cell_labels, "_", as.character(1:length(cell_labels)))
  
  
  #-------------------
  # Simulate true data
  #-------------------
  
  # Generate cell-type specific peaks
  lapply(which_celltypes, function(celltype){
    
    # Apply different rates per cell depending on group label for generating cell-type specific peaks
    n_cells_this_celltype <- sum(cell_labels == celltype)
    counts_celltype <- bulk[,celltype]
    
    # Define probabilities
    #                        Prob observting frag                Total number of fragments epxpected; the 0.5s are for two alleles that will be simulated/added later
    prob_per_peaks <- counts_celltype/sum(counts_celltype) * (n_frags_per_cell*0.5 * (1-rate_noise)) + ((rate_noise*n_frags_per_cell)/n_peaks*0.5) 
    
    # Cap probabilities at something sensible
    prob_per_peaks <- ifelse(prob_per_peaks > 0.9, 0.9, prob_per_peaks)
    
    # Represent the two haplotypes as two random draws
    mat1 <- (matrix(rbinom(n_peaks*n_cells_this_celltype, size = 1, prob = prob_per_peaks),
                    ncol = n_cells_this_celltype, byrow = FALSE) )
    mat2 <- (matrix(rbinom(n_peaks*n_cells_this_celltype, size = 1, prob = prob_per_peaks),
                    ncol = n_cells_this_celltype, byrow = FALSE) )
    
    mat <- mat1 + mat2
    Matrix(mat)
  }) %>% do.call(what = "cbind") -> sparse_matrix
  
  colnames(sparse_matrix) <- final_names
  sparse_matrix
}

# Here, we call the function above to simulate data
simulated_noisy <- simulate_scatac(20, colnames(bulk), rate_noise = 0.2)
simulated_clean <- simulate_scatac(20, colnames(bulk), rate_noise = 0)

# Do a basic LSI embedding to assess
compute_LSI <- function(x){
  nfreqs <- t(t(x) / Matrix::colSums(x))
  idf <- as(log(1 + ncol(x) / Matrix::rowSums(x)), "sparseVector")
  tf_idf_counts <- as(Diagonal(x=as.vector(idf)), "sparseMatrix") %*% nfreqs
  SVD_x <-  irlba(tf_idf_counts, 3, 3)
  d_diag = matrix(0, nrow=length(SVD_x$d), ncol=length(SVD_x$d))
  diag(d_diag) = SVD_x$d
  LSI_x_final = t(d_diag %*% t(SVD_x$v))
  LSI_x_final
}

# Function to do LSI and then create the corresponding data frame
makeLSI_df <- function(simulated){
  # Compute LSI and extract cell types from previous simulation
  LSI_dims <- compute_LSI(simulated)
  celltypes <- str_split_fixed(colnames(simulated), "_", 2)[,1]
  
  # Make one data frame for plotting
  LSI_df <- data.frame(
    LSI_2 = LSI_dims[,2],
    LSI_3 = LSI_dims[,3],
    celltype = celltypes,
    cell_id = colnames(simulated)
  )
  LSI_df
}

# Create two LSI dfs to compare
LSI_df_noise <- makeLSI_df(simulated_noisy)
LSI_df_clean <- makeLSI_df(simulated_clean)

p1 <- ggplot(shuf(LSI_df_clean), aes(x = LSI_2, y = LSI_3, color = celltype)) +
  geom_point(size = 1) + scale_color_manual(values = all_color_maps) +
  ggtitle("clean - simulated")

p2 <- ggplot(shuf(LSI_df_noise), aes(x = LSI_2, y = LSI_3, color = celltype)) +
  geom_point(size = 1) + scale_color_manual(values = all_color_maps) +
  ggtitle("noisy - simulated")

cowplot::ggsave(cowplot::plot_grid(p1, p2, nrow = 1), 
                filename = "../output/erythroid_simulated_comparison.pdf", width = 9, height = 4)




