#!/usr/bin/env Rscript

# Load required libraries
print('Loading required libraries...')
library(SingleCellExperiment)
library(dreamlet)
library(BiocParallel)
library(dplyr)
library(purrr)

# Specify settings
num_cores <- as.integer(Sys.getenv("NSLOTS", "1"))
dir_path <- "/wynton/group/pollen/reedmcmullen/projects/NEMP25/dreamlet_NEMP25/"
data_path <- "/wynton/group/pollen/reedmcmullen/projects/NEMP25/dreamlet_NEMP25/objects/"
res_path <- "/wynton/group/pollen/reedmcmullen/projects/NEMP25/dreamlet_NEMP25/dreamlet_results/"

# Load the dataset
sce = readRDS(paste0(data_path, "NEMP25_Tel_v5_dreamlet.rds"))

# Update SCE object to include group and sample metadata
print('Formatting SCE object...')

levels <- sort(unique(sce$cell_class))
sce$cell_class <- factor(sce$cell_class, levels=levels)
print('Cell class levels set')

levels <- c('Human', 'Chimp', 'Orangutan')
sce$species <- factor(sce$species, levels=levels)
print('Species levels set')

levels <- sort(unique(sce$seq_batch))
sce$seq_batch <- factor(sce$seq_batch, levels=levels)
print('Seq. batch levels set')

levels <- c('6hr', '54hr', '7day')
sce$time_point <- factor(sce$time_point, levels=levels)
print('Time point levels set')

levels <- sort(unique(sce$perturbation))
sce$perturbation <- factor(sce$perturbation, levels=levels)
print('Perturbation levels set')

levels <- sort(unique(sce$individual))
sce$individual <- factor(sce$individual, levels=levels)
print('Individuals levels set')

sce$group <- paste(sce$species, sce$perturbation, sce$time_point, sep='_')
levels <- sort(unique(sce$group))
sce$group <- factor(sce$group, levels=levels)
print('Group levels set')

sce$sample <- paste(sce$group, sce$individual, sep='_')
levels <- sort(unique(sce$sample))
sce$sample <- factor(sce$sample, levels=levels)
print('Sample levels set')

# Define vectors to iterate through for defining dreamlet contrasts
species <- c('Human', 'Chimp')
times <- c('6hr', '54hr', '7day')
treman_perts <- c('FGF8b', 'EGF', 'BMP7', 'IGF1LR3', 'FGF2', 'NRG1', 'SHH', 'JAG1', 'DLL1')
dmso_perts <- c('LY411575', 'HX531', 'LDN193189', 'Rapamycin', 'ATRA', 'BGJ398', 'TRULI', 'WNTC59', 'AZD8931', 'CHIR99021', 'SANT1')
perts <- c(treman_perts, dmso_perts)
cntrls <- c('TreMan', 'DMSO')
baselines <- c('Untreated', 'TreMan', 'DMSO')

# Define response contrasts
responsediff_cons <- c() # Initialize an empty vector
for (spec in species){
  for (time in times){
    for (pert in cntrls){
      con_name <- paste(pert, time, spec, "ResponseDiff", sep = "_")
      con <- paste0("(group", spec, "_", pert, "_", time, " - group", spec, "_Untreated_", time, ")")
      responsediff_cons[con_name] <- con
    }
    for (pert in treman_perts){
      con_name <- paste(pert, time, spec, "ResponseDiff", sep = "_")
      con <- paste0("(group", spec, "_", pert, "_", time, " - group", spec, "_TreMan_", time, ")")
      responsediff_cons[con_name] <- con
    }
    for (pert in dmso_perts){
      con_name <- paste(pert, time, spec, "ResponseDiff", sep = "_")
      con <- paste0("(group", spec, "_", pert, "_", time, " - group", spec, "_DMSO_", time, ")")
      responsediff_cons[con_name] <- con
    }
  }
}

print(paste0("N response difference contrasts: ", length(responsediff_cons)))

# Downsampling analysis
n_iterations <- 500

# Get cell metadata
cell_metadata <- as.data.frame(colData(sce))

# Generate an empty list to store results
all_iteration_results <- list()

for (iteration in seq_len(n_iterations)) {
  message("Running iteration ", iteration, "/", n_iterations)

  # Reproducible random sampling for this iteration
  set.seed(41 + iteration)  
  keep_cells <- character(0)



  # Loop over species
  for (species_group in unique(cell_metadata$group)) {
    print('Downsampling...')
    # Cells belonging to this group
    group_cells <- rownames(cell_metadata)[cell_metadata$group == species_group]

    # Human and Chimp cells within the group
    chimp <- group_cells[cell_metadata[group_cells, "species"] == "Chimp"]
    human <- group_cells[cell_metadata[group_cells, "species"] == "Human"]

    # Match to the smaller species count
    n <- min(length(chimp), length(human))

    # Only sample if both species are present
    if (n > 0) {
      keep_cells <- c(
        keep_cells,
        sample(chimp, size = n, replace = FALSE),
        sample(human, size = n, replace = FALSE)
      )
    }
  }

  # Create balanced SCE for this iteration
  sce_balanced <- sce[, keep_cells]
  # Pseudobulk data based on sample
  print('Pseudobulking...')
  pb <- aggregateToPseudoBulk(sce_balanced,
                            assay="X",
                            cluster_id="cell_class",
                            sample_id="sample",
                            verbose=TRUE)
  # Normalize and apply voomWithDreamWeights
  print('Running voomWithDreamWeights...')
  model_form = ~ (1|species) + (1|perturbation) + (1|time_point) + (1|individual)
  res_voom = processAssays(pb,
                         model_form,
                         min.count=5, # Drop genes with fewer pseudobulk counts
                         min.cells=50, # Drop samples with fewer cells
                         min.samples=2, # Minimum number of samples passing cutoffs for cell cluster to be retained
                         min.prop=0.00483, # Minimum proportion of retained samples with non-zero counts for a gene to be retained
                         BPPARAM=SnowParam(num_cores, type = "SOCK"),
                         quiet=FALSE)
  # Run dreamlet
  print('Running dreamlet...')
  dreamlet_form = ~ 0 + group
  res_dreamlet = dreamlet(res_voom,
                  dreamlet_form,
                  contrasts = responsediff_cons,
                  BPPARAM = SnowParam(num_cores, type="SOCK", progressbar=TRUE))

  # Format dreamlet results
  print('Formatting dreamlet results...')
  contrasts <- coefNames(res_dreamlet)
  results_df = contrasts %>%
  purrr::map_df(function(x)
    topTable(res_dreamlet, coef=x, number=Inf) %>%
    as_tibble %>%
    mutate(coef=factor(x, contrasts))) %>%
    mutate(se = logFC / t) %>% # Compute standard error and add to se column.
    mutate(significant = ifelse(adj.P.Val < 0.10, TRUE, FALSE)) %>% # Add significance flag.
    mutate(logFCsign = sign(logFC)) # Add logFCsign column with sign of logFC.
  
  # Store results
  all_iteration_results[[iteration]] <- results_df
}

# Combine results into one data frame
print('Combining results...')
all_results <- dplyr::bind_rows(all_iteration_results)

# Save results
saveRDS(all_results, paste0(res_path, 'dreamlet_response_downsampling_', n_iterations, 'iterations_results_df.rds'))
print('Finished.')



