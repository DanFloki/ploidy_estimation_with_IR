# AUTHOR: Daniele Buono
# daniele.buono@tu-dresden.de
# daniele.buo@gmail.com
# last edit: 31 May 2022

# How to run this script #####
# The script does not run automatically as a whole, but each line has to be executed individually (recommended R-studio)
# reading and filtering raw data has been performed separately for herbarium (dry) 
# and living (fresh) specimens, using the same script (this) but changing variable and folder names.
# Combination of the two dataset (combined) has been performed afterwards.
# Follow filtering steps from 0 to 6:
#   0. read spectra
#   1. remove spectra that have values > 100 % reflectance
#   2. remove ID_sample with less than 3 spectra
#   3. apply Savitzky - Golay filter
#   4. filter spectra too similar to background
#   5. separate samples with unknown ploidy from samples with known ploidy
#   6. reduce dataset size by averaging value every 5 nm
# Afterwards, some elaboration on the datasets:
#   Average by sample ID
#   Subset by sampling point
#   Average by individual ID
#   combining herbarium and fresh dataset
# Then to the analysis:
#   plot spectra
#   hiearchical cluster analysis
#   DAPC
#   NeuralNetwork - NOTE: optimization of parameters for NNET performed on HPC cluster using other scripts
#   Statistics on DAPC and NeuralNetwork results
#
# Before running the instructions it is necessary to create the folder structure indicated in the various setwd and in lines 24-27 of IR_spectroscopy_functions.R


# set working directory
setwd("~/Dropbox/Infrared Spectroscopy/R script analysis")

# load sources #####
# NOTE: before continue, open IR_spectroscopy_functions.R and change spectra_dataset value in:
# spectra_dataset <-"dry" # to use herbarium spectra
# spectra_dataset <- "fresh" # to use spectra of fresh (living) plants
# spectra_dataset <- "combined_dry_fresh" # to use spectra of herbarium and live


source("IR_spectroscopy_functions.R")

#read reference table ####
# non filtered data, load one of the two:
raw_ref_table <- read_reference_table("0_raw_herbarium_reference_table.csv", dir="filt") # see note above
#raw_ref_table <- read_reference_table("0_fresh_Veronica_refdata.csv", dir="filt") # see note above


summary(raw_ref_table) #general summary on dataset
unique(raw_ref_table$ID_ind) #number of unique specimens (individuals)
nrow(raw_ref_table)
head(raw_ref_table)

# 0. read the spectra from several output files and create a single matrix to store them ####
# raw files availabe at https://zenodo.org/records/7729563
setwd("dry/spectra dry material/lon_spi_herbarium/raw_data") # for herbarium material
# setwd("fresh/raw_data") #for fresh material

raw_spectra_matrix <- read_raw_spectrophotometer_data(raw_ref_table)
e <- edit(raw_spectra_matrix)
#save spectra as matrix on file
setwd("/home/daniele/Dropbox/Infrared Spectroscopy/R script analysis/dry/filtered_data/")
# setwd("/home/daniele/Dropbox/Infrared Spectroscopy/R script analysis/fresh/filtered_data/")

write.table(raw_spectra_matrix,"0_raw_herbarium_spectra.txt", sep="\t")
# write.table(raw_spectra_matrix,"0_raw_fresh_spectra.txt",sep="\t")

# 1. remove spectra that have values > 100 % reflectance ####
raw_spectra_matrix<-read_spect_matrix("0_raw_herbarium_spectra.txt", dir="filt") # for herbarium material
#raw_spectra_matrix<-read_spect_matrix("0_raw_fresh_spectra.txt", dir="filt") # for fresh material
filt_1_spect_matrix <- remove_anomalous_spectra(raw_spectra_matrix)
filt_1_ref_table <- raw_ref_table[raw_ref_table$ID_spect %in% rownames(filt_1_spect_matrix), ]
#check if reference table corresponds to spectra in matrix
identical(as.character(filt_1_ref_table$ID_spect),row.names(filt_1_spect_matrix))

# herbarium:
write.table(filt_1_ref_table,"1_ref_table_herbarium_without_abnormal.txt",sep="\t")
write.table(filt_1_spect_matrix,"1_spectra_matrix_herbarium_without_abnormal.txt",sep="\t")
# fresh:
# write.table(filt_1_ref_table,"1_ref_table_fresh_without_abnormal.txt",sep="\t") 
# write.table(filt_1_spect_matrix,"1_spectra_matrix_fresh_without_abnormal.txt",sep="\t")

# 2. remove ID_sample with less than 3 spectra ####
# herbarium:
filt_1_spect_matrix<-read_spect_matrix("1_spectra_matrix_herbarium_without_abnormal.txt", dir="filt")
filt_1_ref_table <- read_reference_table("1_ref_table_herbarium_without_abnormal.txt", dir="filt")
# fresh:
# filt_1_spect_matrix<-read_spect_matrix("1_spectra_matrix_fresh_without_abnormal.txt", dir="filt")
# filt_1_ref_table <- read_reference_table("1_ref_table_fresh_without_abnormal.txt", dir="filt")

filt_2_spect_matrix <- remove_sampled_less_than_3(filt_1_ref_table, filt_1_spect_matrix)
filt_2_ref_table <- filt_1_ref_table[filt_1_ref_table$ID_spect %in% rownames(filt_2_spect_matrix), ]
nrow(filt_1_ref_table)
nrow(filt_2_ref_table)
nrow(filt_2_spect_matrix)
identical(as.character(filt_2_ref_table$ID_spect),row.names(filt_2_spect_matrix)) # check if TRUE
write.table(filt_2_spect_matrix, "2_spectra_matrix_more_than_2_ID_sample.txt", sep="\t")
write.table(filt_2_ref_table, "2_ref_table_more_2_ID_sample.txt", sep="\t")

# 3. Savitzky - Golay filter ####
filt_2_spect_matrix <- read_spect_matrix("2_spectra_matrix_more_than_2_ID_sample.txt", dir="filt")
filt_2_ref_table <- read_reference_table("2_ref_table_more_2_ID_sample.txt", dir="filt")

filt_3_spect_matrix <- apply_sgolay_filter(filt_2_spect_matrix, 
                                           p=3, n=33) # p=3 and n=33 values assigned by looking at graph of real spectra (see function)

write.table(filt_3_spect_matrix, "3_spectra_matrix_SGolay.txt", sep="\t")
write.table(filt_2_ref_table, "3_filt_ref_table_SGolay.txt", sep="\t")

# 4. filter spectra too similar to background ####
# NOTE: this filtering step only for herbarium material
filt_3_spect_matrix <- read_spect_matrix("3_spectra_matrix_SGolay.txt", dir="filt")
filt_3_ref_table <- read_reference_table("3_filt_ref_table_SGolay.txt", dir="filt")

background <- as.matrix(read.table("4_avg_spect_matrix_S-Golay_paper.txt", sep="\t", header = T))

filt_4_spect_matrix <- remove_similar_to_background(spect_matrix = filt_3_spect_matrix, 
                             ref_table = filt_3_ref_table, 
                             spect_matrix_background = background, 
                             k=10)
filt_4_ref_table <- filt_3_ref_table[filt_3_ref_table$ID_spect %in% rownames(filt_4_spect_matrix), ]

nrow(filt_4_spect_matrix)
nrow(filt_4_ref_table)
identical(as.character(filt_4_ref_table$ID_spect),row.names(filt_4_spect_matrix))

write.table(filt_4_spect_matrix,"4_spect_matrix_SGolay_without_similar_background_k10.txt", sep="\t")
write.table(filt_4_ref_table,"4_ref_table_SGolay_without_similar_background_k10.txt", sep="\t")

# 5. separate samples with unknown ploidy from samples with known ploidy (used to create statistic model) ####
# NOTE: for living individuals we know the ploidy levels of all individuals
filt_4_spect_matrix<-read_spect_matrix("4_spect_matrix_SGolay_without_similar_background_k10.txt", dir="filt")
filt_4_ref_table <- read_reference_table("4_ref_table_SGolay_without_similar_background_k10.txt", dir="filt")

filt_5_spect_matrix <- filt_4_spect_matrix[!is.na(filt_4_ref_table$ploidy),]
filt_5_ref_table<-filt_4_ref_table[!is.na(filt_4_ref_table$ploidy),]
identical(as.character(filt_5_ref_table$ID_spect),row.names(filt_5_spect_matrix))

nrow(filt_5_spect_matrix)
nrow(filt_5_ref_table)

write.table(filt_5_spect_matrix,"5_spect_matrix_SGolay_background_more2sample_ploidy.txt",sep="\t")
write.table(filt_5_ref_table,"5_ref_table_SGolay_background_more2sample_ploidy.txt",sep="\t")

filt_5_spect_matrix_UNK_PLO <- filt_4_spect_matrix[is.na(filt_4_ref_table$ploidy),]
filt_5_ref_table_UNK_PLO<-filt_4_ref_table[is.na(filt_4_ref_table$ploidy),]
write.table(filt_5_spect_matrix_UNK_PLO,"5_UNKNOWN_PLOIDY_spect_matrix_SGolay_background_more2sample.txt",sep="\t")
write.table(filt_5_ref_table_UNK_PLO,"5_UNKNOWN_PLOIDY_ref_table_SGolay_background_more2sample.txt",sep="\t")

# 6. reduce dataset size by averaging value every 5 nm ###### 
filt_5_spect_matrix <- read_spect_matrix("5_spect_matrix_SGolay_background_more2sample_ploidy.txt", dir="filt")
filt_5_ref_table <- read_reference_table("5_ref_table_SGolay_background_more2sample_ploidy.txt", dir="filt")

# NOTE: for fresh, we skip filtering steps 4 and 5
filt_5_spect_matrix <- read_spect_matrix("3_spectra_matrix_SGolay.txt", dir="filt")
filt_5_ref_table <- read_reference_table("3_filt_ref_table_SGolay.txt", dir="filt")


filt_6_spect_matrix <- reduce_dataset_size(spect_matrix = filt_5_spect_matrix, 
                                           reduction_factor = 5)

e <- edit(filt_6_spect_matrix)

write.table(filt_6_spect_matrix, "6_spect_matrix_reduced_dataset.txt",sep="\t")

# filt_4: spectra filtered
setwd("~/Dropbox/Infrared Spectroscopy/R script analysis/fresh/filtered_data/filt_3_spectra_plotted/")
col_graphs="#D95F02"
plot_spectra(ref_table = filt_3_ref_table, 
             spect_matrix = filt_3_spect_matrix, 
             col_graphs = col_graphs, 
             max_y = 100)

# filt_6: spectra filtered and reduced in size (avg every 5 nm)
setwd("~/Dropbox/Infrared Spectroscopy/R script analysis/fresh/filtered_data/filt_6_spectra_plotted/")
plot_spectra(ref_table = filt_3_ref_table, # ref_tage for filt_6 dataset is the same as for filt_5 dataset
             spect_matrix = filt_6_spect_matrix, 
             col_graphs = col_graphs, 
             max_y = 100)

# count specimens with less than 24 sampling points ####
filt_ref_table <- read_reference_table("ref_table_filt_6.txt", dir="filt")
filt_spect_matrix <- read_spect_matrix("spectra_matrix_filt_6.txt", dir="filt")
unique(filt_ref_table$ID_ind)


sampled_24=0
sampled_less_24=0
sampled_at_least_12=0
sampled_at_least_20=0
sampled_at_least_6=0
for(i in unique(filt_ref_table$ID_ind)){
  # specimens with 24 sampling points
  if(nrow(filt_ref_table[filt_ref_table$ID_ind == i,])==24){
    sampled_24=sampled_24+1
  }
  
  #specimens with at least 12 sampling points
  if (nrow(filt_ref_table[filt_ref_table$ID_ind == i,])>11){
    sampled_at_least_12=sampled_at_least_12+1
  }
  
  #specimens with at least 6 sampling points
  if (nrow(filt_ref_table[filt_ref_table$ID_ind == i,])>5){
    sampled_at_least_6=sampled_at_least_6+1
  }
  
  #specimens with at least 6 sampling points
  if (nrow(filt_ref_table[filt_ref_table$ID_ind == i,])>19){
    sampled_at_least_20=sampled_at_least_20+1
  }
  
  #specimens with less than 24 sampling points
  if (nrow(filt_ref_table[filt_ref_table$ID_ind == i,])<24){
    sampled_less_24=sampled_less_24+1
  }
}

n_specimens = nlevels(filt_ref_table$ID_ind)
n_specimens
sampled_24
sampled_at_least_12
sampled_at_least_6
sampled_less_24
sampled_at_least_20

# average by ID_sample ####
filt_6_spect_matrix <- read_spect_matrix("6_spect_matrix_reduced_dataset.txt", dir="filt")
#filt_6_ref_table <- read_reference_table("5_ref_table_SGolay_background_more2sample_ploidy.txt", dir="filt")
filt_6_ref_table <- read_reference_table("3_filt_ref_table_SGolay.txt", dir="filt")


sample_avg_spectra_matrix <- average_sample_spectra(filt_6_ref_table, 
                                                    filt_6_spect_matrix)

nrow(sample_avg_spectra_matrix$ref_table)
nrow(sample_avg_spectra_matrix$spect_matrix)

setwd("~/Dropbox/Infrared Spectroscopy/R script analysis/fresh/elaborated_data/plots_average_spectra_ID_sample_filt_6/")
plot_ind_spectra(spect_matrix = sample_avg_spectra_matrix$spect_matrix,
                 col_graphs = col_graphs, max_y = 100)

setwd("~/Dropbox/Infrared Spectroscopy/R script analysis/fresh/elaborated_data/")

write.table(sample_avg_spectra_matrix$spect_matrix, "average_spectra_ID_sample_filt_6.txt",sep="\t")
write.table(sample_avg_spectra_matrix$ref_table, "ref_table_average_spectra_ID_sample_filt_6.txt", sep="\t")

# subset dataset by sampling point ####
filt_6_spect_matrix <- read_spect_matrix("spectra_matrix_filt_6.txt", dir="filt")
filt_6_ref_table <- read_reference_table("ref_table_filt_6.txt", dir="filt")
nrow(filt_6_spect_matrix)
nrow(filt_6_ref_table)
identical(as.character(filt_6_ref_table$ID_spect),row.names(filt_6_spect_matrix))

nlevels(filt_6_ref_table$ID_ind)

filt_6_spect_matrix_bottom <- filt_6_spect_matrix[filt_6_ref_table$leaf_side == "b", ]
filt_6_ref_table_bottom <- filt_6_ref_table[filt_6_ref_table$leaf_side == "b", ]


identical(as.character(filt_6_ref_table_bottom$ID_spect),row.names(filt_6_spect_matrix_bottom))
nrow(filt_6_spect_matrix_bottom)
nrow(filt_6_ref_table_bottom)


write.table(filt_6_spect_matrix_bottom,"spectra_matrix_filt_6_bottom.txt",sep="\t")
write.table(filt_6_ref_table_bottom,"ref_table_filt_6_bottom.txt",sep="\t")


filt_6_spect_matrix_top <- filt_6_spect_matrix[filt_6_ref_table$leaf_side == "t", ]
filt_6_ref_table_top <- filt_6_ref_table[filt_6_ref_table$leaf_side == "t", ]

identical(as.character(filt_6_ref_table_top$ID_spect),row.names(filt_6_spect_matrix_top))
nrow(filt_6_spect_matrix_top)
nrow(filt_6_ref_table_top)

write.table(filt_6_spect_matrix_top,"spectra_matrix_filt_6_top.txt",sep="\t")
write.table(filt_6_ref_table_top,"ref_table_filt_6_top.txt",sep="\t")

# subset on avg ID sample
filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_sample_filt_6.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_sample_filt_6.txt", dir="el")
identical(as.character(filt_ref_table$ID_spect),row.names(filt_spect_matrix))
nrow(filt_ref_table)

filt_6_spect_matrix_ID_sample_bottom <- filt_spect_matrix[filt_ref_table$leaf_side == "b", ]
filt_6_ref_table_ID_sample_bottom <- filt_ref_table[filt_ref_table$leaf_side == "b", ]

identical(as.character(filt_6_ref_table_ID_sample_bottom$ID_spect),row.names(filt_6_spect_matrix_ID_sample_bottom))
nrow(filt_6_ref_table_ID_sample_bottom)

write.table(filt_6_spect_matrix_ID_sample_bottom,"spectra_matrix_filt_6_ID_sample_bottom.txt",sep="\t")
write.table(filt_6_ref_table_ID_sample_bottom,"ref_table_filt_6_ID_sample_bottom.txt",sep="\t")


filt_6_spect_matrix_ID_sample_top <- filt_spect_matrix[filt_ref_table$leaf_side == "t", ]
filt_6_ref_table_ID_sample_top <- filt_ref_table[filt_ref_table$leaf_side == "t", ]

identical(as.character(filt_6_ref_table_ID_sample_top$ID_spect),row.names(filt_6_spect_matrix_ID_sample_top))

write.table(filt_6_spect_matrix_ID_sample_top,"spectra_matrix_filt_6_ID_sample_top.txt",sep="\t")
write.table(filt_6_ref_table_ID_sample_top,"ref_table_filt_6_ID_sample_top.txt",sep="\t")

# subset on avg ID Individual
filt_6_spect_matrix_bottom <- read_spect_matrix("spectra_matrix_filt_6_bottom.txt", dir="filt")
filt_6_ref_table_bottom <- read_reference_table("ref_table_filt_6_bottom.txt", dir="filt")

indiv_avg_spectra_matrix_bottom <- average_individual_spectra(ref_table = filt_6_ref_table_bottom, 
                                                       spect_matrix = filt_6_spect_matrix_bottom)

identical(as.character(indiv_avg_spectra_matrix_bottom$ref_table$ID_spect),
          row.names(indiv_avg_spectra_matrix_bottom$spect_matrix))
nrow(indiv_avg_spectra_matrix_bottom$ref_table)
nrow(indiv_avg_spectra_matrix_bottom$spect_matrix)

col_graphs="#D95F02"
plot_ind_spectra(spect_matrix = indiv_avg_spectra_matrix_bottom$spect_matrix,
                 col_graphs = col_graphs, max_y = 100)

write.table(indiv_avg_spectra_matrix_bottom$spect_matrix, "average_spectra_ID_individual_filt_6_bottom.txt", sep="\t")
write.table(indiv_avg_spectra_matrix_bottom$ref_table, "ref_table_average_spectra_ID_individual_filt_6_bottom.txt", sep="\t")



filt_6_spect_matrix_top <- read_spect_matrix("spectra_matrix_filt_6_top.txt", dir="filt")
filt_6_ref_table_top <- read_reference_table("ref_table_filt_6_top.txt", dir="filt")

indiv_avg_spectra_matrix_top <- average_individual_spectra(ref_table = filt_6_ref_table_top, 
                                                              spect_matrix = filt_6_spect_matrix_top)

identical(as.character(indiv_avg_spectra_matrix_top$ref_table$ID_spect),row.names(indiv_avg_spectra_matrix_top$spect_matrix))


col_graphs="#D95F02"
plot_ind_spectra(spect_matrix = indiv_avg_spectra_matrix_top$spect_matrix,
                 col_graphs = col_graphs, max_y = 100)

write.table(indiv_avg_spectra_matrix_top$spect_matrix, "average_spectra_ID_individual_filt_6_top.txt", sep="\t")
write.table(indiv_avg_spectra_matrix_top$ref_table, "ref_table_average_spectra_ID_individual_filt_6_top.txt", sep="\t")



# average by ID individual ####
filt_6_spect_matrix <- read_spect_matrix("spectra_matrix_filt_6.txt", dir="filt")
filt_6_ref_table <- read_reference_table("ref_table_filt_6.txt", dir="filt")

indiv_avg_spectra_matrix <- average_individual_spectra(ref_table = filt_6_ref_table, 
                                                       spect_matrix = filt_6_spect_matrix)

col_graphs="#D95F02"
setwd("~/Dropbox/Infrared Spectroscopy/R script analysis/fresh/elaborated_data/plots_average_spectra_ID_ind_filt_6/")
plot_ind_spectra(spect_matrix = indiv_avg_spectra_matrix$spect_matrix,
                 col_graphs = col_graphs, max_y = 100)

setwd("~/Dropbox/Infrared Spectroscopy/R script analysis/fresh/elaborated_data/")

write.table(indiv_avg_spectra_matrix$spect_matrix, "average_spectra_ID_individual_filt_6.txt", sep="\t")
write.table(indiv_avg_spectra_matrix$ref_table, "ref_table_average_spectra_ID_individual_filt_6.txt", sep="\t")

# for dataset subsetted by sampling point
indiv_avg_spectra_matrix_bottom <- average_individual_spectra(ref_table = filt_6_ref_table_bottom, 
                                                       spect_matrix = filt_6_spect_matrix_bottom)

write.table(indiv_avg_spectra_matrix_bottom$spect_matrix, "average_spectra_ID_individual_filt_6_bottom.txt", sep="\t")
write.table(indiv_avg_spectra_matrix_bottom$ref_table, "ref_table_average_spectra_ID_individual_filt_6_bottom.txt", sep="\t")

indiv_avg_spectra_matrix_top <- average_individual_spectra(ref_table = filt_6_ref_table_top, 
                                                              spect_matrix = filt_6_spect_matrix_top)

write.table(indiv_avg_spectra_matrix_top$spect_matrix, "average_spectra_ID_individual_filt_6_top.txt", sep="\t")
write.table(indiv_avg_spectra_matrix_top$ref_table, "ref_table_average_spectra_ID_individual_filt_6_top.txt", sep="\t")


# combining herbarium and fresh datasets ####
# 1. set up spectra_dataset <-"dry" in IR_spectroscopy_functions.R
dry_ref_table <- read_reference_table("5_ref_table_SGolay_background_more2sample_ploidy.txt", dir="filt")
dry_spect_matrix <- read_spect_matrix("6_spect_matrix_reduced_dataset.txt", dir="filt")
# 2. set up spectra_dataset <-"fresh" in IR_spectroscopy_functions.R
fresh_ref_table <- read_reference_table("6_filt_ref_table_reduced_dataset.txt", dir="filt")
fresh_spect_matrix <- read_spect_matrix("6_spect_matrix_reduced_dataset.txt", dir="filt")

#3. combining the datasets (NOTE: using rbind, does not matter if the order of variables in the dataset is different)
# (first removing a column that is only in dry dataset)
dry_ref_table_2<-dry_ref_table[,-1]
dry_ref_table_2
comb_ref_table <- rbind(dry_ref_table_2, fresh_ref_table)

fresh_ref_table
View(comb_ref_table)
comb_spect_matrix <- rbind(dry_spect_matrix, fresh_spect_matrix)

#4. setting the working directory and save joined datasets
setwd("~/Dropbox/Infrared Spectroscopy/R script analysis/combined_dry_fresh/filtered_data")

write.table(comb_spect_matrix, "spectra_matrix_filt_6.txt", sep="\t")
write.table(comb_ref_table, "ref_table_filt_6.txt", sep="\t")

#same for average_ID_ind
# 1. set up spectra_dataset <-"dry" in IR_spectroscopy_functions.R
dry_ref_table <- read_reference_table("ref_table_average_spectra_ID_individual_filt_6.txt", dir="el")
dry_spect_matrix <- read_spect_matrix("average_spectra_ID_individual_filt_6.txt", dir="el")
# 2. set up spectra_dataset <-"fresh" in IR_spectroscopy_functions.R
fresh_ref_table <- read_reference_table("ref_table_average_spectra_ID_individual_filt_6.txt", dir="el")
fresh_spect_matrix <- read_spect_matrix("average_spectra_ID_individual_filt_6.txt", dir="el")

#3. combining the datasets (NOTE: using rbind, does not matter if the order of variables in the dataset is different)

comb_ref_table <- rbind(dry_ref_table, fresh_ref_table)
View(comb_ref_table)
comb_spect_matrix <- rbind(dry_spect_matrix, fresh_spect_matrix)

#4. setting the working directory and save joined datasets
setwd("~/Dropbox/Infrared Spectroscopy/R script analysis/combined_dry_fresh/elaborated_data")

write.table(comb_spect_matrix, "average_spectra_ID_individual_filt_6.txt", sep="\t")
write.table(comb_ref_table, "ref_table_average_spectra_ID_individual_filt_6.txt", sep="\t")

# plot example spectra ####
# plotting herbarium OLD01319 and living lo32aUkr_t: they should be the same individual
# ALL spectra
filt_ref_table <- read_reference_table("ref_table_filt_6.txt", dir="filt")
filt_spect_matrix <- read_spect_matrix("spectra_matrix_filt_6.txt", dir="filt")


two_indiv_ref <- filt_ref_table[filt_ref_table$ID_ind == "OLD01319" | filt_ref_table$ID_ind == "lo32aUkr_t",]
two_indiv_spect <- filt_spect_matrix[filt_ref_table$ID_ind == "OLD01319" | filt_ref_table$ID_ind == "lo32aUkr_t",]
setwd("~/Dropbox/Infrared Spectroscopy/manuscript/figures/for scheme/")
col_graphs="#D95F02"

plot_spectra_std_err(ref_table = two_indiv_ref, 
             spect_matrix = two_indiv_spect, 
             col_graphs = col_graphs, 
             max_y = 100)

plot_avg_spectra_std_err(ref_table = two_indiv_ref, 
                     spect_matrix = two_indiv_spect, 
                     col_graphs = col_graphs, 
                     max_y = 100)

plot_avg_spectra_std_err_one_graph(ref_table = two_indiv_ref, 
                         spect_matrix = two_indiv_spect,
                         max_y = 100,
                         file_name = "two individuals spectra in one plot")
# plot average 2x and 4x separately ####
# for fresh, change line 22 in IR_spectroscopy_funcitons.R:  spectra_dataset <- "fresh"
filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_individual_filt_6.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_individual_filt_6.txt", dir="el")

spicata_spectra <- filt_spect_matrix[filt_ref_table$species == "spicata",]
spicata_ref_table <- filt_ref_table[filt_ref_table$species == "spicata",]

longifolia_spectra <- filt_spect_matrix[filt_ref_table$species == "longifolia",]
longifolia_ref_table <- filt_ref_table[filt_ref_table$species == "longifolia",]



setwd("~/Dropbox/Infrared Spectroscopy/R script analysis/plots")

col_graphs=c("#D95F02", "#1B9E77")

plot_avg_and_std_err_separation_ploidy(spicata_ref_table, spicata_spectra, col_graphs, 100,"spicata_fresh.svg")
plot_avg_and_std_err_separation_ploidy(longifolia_ref_table, longifolia_spectra, col_graphs, 100,"longifolia_fresh.svg")

# for dry samples, uncomment line 20 in IR_spectroscopy_functions.R!

filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_individual_filt_6.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_individual_filt_6.txt", dir="el")

spicata_spectra <- filt_spect_matrix[filt_ref_table$species == "spicata",]
spicata_ref_table <- filt_ref_table[filt_ref_table$species == "spicata",]

longifolia_spectra <- filt_spect_matrix[filt_ref_table$species == "longifolia",]
longifolia_ref_table <- filt_ref_table[filt_ref_table$species == "longifolia",]



setwd("~/Dropbox/Infrared Spectroscopy/R script analysis/plots")

col_graphs=c("#D95F02", "#1B9E77")

plot_avg_and_std_err_separation_ploidy(spicata_ref_table, spicata_spectra, col_graphs, 100,"spicata_dry.svg")
plot_avg_and_std_err_separation_ploidy(longifolia_ref_table, longifolia_spectra, col_graphs, 100,"longifolia_dry.svg")

# plots used for scheme
one_ind_ref<-two_indiv_ref[two_indiv_ref$ID_ind == "lo32aUkr_t",]
one_ind_spect<-two_indiv_spect[two_indiv_ref$ID_ind == "lo32aUkr_t",]

top_one_ind_ref<-one_ind_ref[one_ind_ref$leaf_side=="t",]
top_one_ind_spect<-one_ind_spect[one_ind_ref$leaf_side=="t",]
bottom_one_ind_ref<-one_ind_ref[one_ind_ref$leaf_side=="b",]
bottom_one_ind_spect<-one_ind_spect[one_ind_ref$leaf_side=="b",]

plot_spectra_std_err(ref_table = top_one_ind_ref, 
                     spect_matrix = top_one_ind_spect, 
                     col_graphs = col_graphs, 
                     max_y = 100,
                     avg_std_err_plot = F)
plot_spectra_std_err(ref_table = bottom_one_ind_ref, 
                     spect_matrix = bottom_one_ind_spect, 
                     col_graphs = col_graphs, 
                     max_y = 100,
                     avg_std_err_plot = F)


plot_avg_spectra_std_err(ref_table = top_one_ind_ref, 
                         spect_matrix = top_one_ind_spect, 
                         col_graphs = col_graphs, 
                         max_y = 100)
plot_avg_spectra_std_err(ref_table = bottom_one_ind_ref, 
                         spect_matrix = bottom_one_ind_spect, 
                         col_graphs = col_graphs, 
                         max_y = 100)

plot_avg_spectra_std_err(ref_table = one_ind_ref, 
                         spect_matrix = one_ind_spect, 
                         col_graphs = col_graphs, 
                         max_y = 100)

# hierarchical cluster analysis ####
filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_sample_filt_6.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_sample_filt_6.txt", dir="el")

filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_individual_filt_6.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_individual_filt_6.txt", dir="el")



dist_matrix = dist(filt_spect_matrix)
h_clust <- hclust(dist_matrix, method = "ward.D2")
cols = c("#1B9E77","#D95F02")

library(dendextend)
phylo_h_clust = as.dendrogram(h_clust)
plot(phylo_h_clust)
phylo_h_clust %>% plot
labels(phylo_h_clust)
# Change color and size for labels
col_labs = character()
filt_ref_table$ploidy <- as.factor(filt_ref_table$ploidy)
for(i in labels(phylo_h_clust)){
  col_i <- cols[filt_ref_table[filt_ref_table$ID_ind==i,]$ploidy]
  col_labs <- c(col_labs, col_i)
}

phylo_h_clust %>% set("labels_col", col_labs) %>% # change color
  set("labels_cex", 0.7) %>% # Change size
  plot(main = "Hierarchical clustering avg ID_individual") # plot
legend(x=80, y=3000, 
       bty="n",# title="ploidy", 
       legend=c(paste0(levels(filt_ref_table$ploidy),"x")), fill=cols)


# DAPC analysis ####

# load one of the following dataset, to choose between herbarium, fresh, or combined, go on line 16 of script IR_spectroscopy_functions.R
# ID sample spectra
filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_sample_filt_6.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_sample_filt_6.txt", dir="el")
filt_ref_table <- read_reference_table("ref_table_filt_6_ID_sample_bottom.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("spectra_matrix_filt_6_ID_sample_bottom.txt", dir="el")
filt_ref_table <- read_reference_table("ref_table_filt_6_ID_sample_top.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("spectra_matrix_filt_6_ID_sample_top.txt", dir="el")


# ID ind spectra (ID SPECIM)
filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_individual_filt_6.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_individual_filt_6.txt", dir="el")
filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_individual_filt_6_bottom.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_individual_filt_6_bottom.txt", dir="el")
filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_individual_filt_6_top.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_individual_filt_6_top.txt", dir="el")

# ALL spectra
filt_ref_table <- read_reference_table("ref_table_filt_6.txt", dir="filt")
filt_spect_matrix <- read_spect_matrix("spectra_matrix_filt_6.txt", dir="filt")
filt_ref_table <- read_reference_table("ref_table_filt_6_bottom.txt", dir="filt")
filt_spect_matrix <- read_spect_matrix("spectra_matrix_filt_6_bottom.txt", dir="filt")
filt_ref_table <- read_reference_table("ref_table_filt_6_top.txt", dir="filt")
filt_spect_matrix <- read_spect_matrix("spectra_matrix_filt_6_top.txt", dir="filt")


# 1. find best number of PCs 
best_PCs_value <- DAPC_find_best_PCs(ref_table = filt_ref_table, spect_matrix = filt_spect_matrix,
                                min_PCA = 20, max_PCA = 20, step_increase = 10,
                                rep_i = nlevels(filt_ref_table$ID_ind) # to test with all individuals
                                )
unique(filt_ref_table$ID_ind)
# perform dapc with best number of PCs found
dapc1 <- dapc(filt_spect_matrix, filt_ref_table$ploidy, n.pca = best_PCs_value, n.da=1)

scatter(dapc1, # col = other(gi.lonspi.filt)$comparePal, 
        cex = 0.5, 
        legend = TRUE, posi.leg = "topright", cleg = 0.8, txt.leg=c("2x","4x"),
        clabel = F,
        scree.pca = T, posi.pca = "topleft", 
        xax = 1, yax = 1, inset.solid = .5)

title(main=paste0("DAPC analysis fresh data, ", best_PCs_value," PCs;
                  dataset = ID sample fresh bottom"))



# perform dapc with best number of PCs found
filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_individual_filt_6.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_individual_filt_6.txt", dir="el")
best_PCs_value = 40
dapc2 <- dapc(filt_spect_matrix, filt_ref_table$ploidy, n.pca = best_PCs_value, n.da=1)

scatter(dapc2, # col = other(gi.lonspi.filt)$comparePal, 
        cex = 0.5, legend = TRUE,
        clabel = FALSE, posi.leg = "bottomleft", scree.pca = F,
        posi.pca = "topleft", cleg = 0.75, xax = 1, yax = 2, inset.solid = 1)
title(main=paste0("DAPC analysis herbarium data, ",best_PCs_value," PCs\n
                  dataset = avg_ID_sample"))

# plot best DAPC results

par(mfrow=c(2,2)) # to have 2 plots in one image
# top 2 best performances:
# 84,1% accuracy, 50 PCs using herbarium specimens, avg by ID individual (ID specimen), top side of leaves
# 70,8% accuracy, 20 PCs using fresh specimens, avg by ID sample, bottom side of leaves

# first, go to IR_spectroscopy_functions.R and change dataset in dry (line 17)
filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_individual_filt_6_top.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_individual_filt_6_top.txt", dir="el")

# ID_speciment, herbarium, top
best_PCs_value = 50
dapc1 <- dapc(filt_spect_matrix, filt_ref_table$ploidy, n.pca = best_PCs_value, n.da=1)

scatter(dapc1, # col = other(gi.lonspi.filt)$comparePal, 
        cex = 0.5, 
        legend = TRUE, posi.leg = "topright", cleg = 0.8, txt.leg=c("2x","4x"),
        clabel = F,
        scree.pca = T, posi.pca = "topleft", 
        xax = 1, yax = 1, inset.solid = .5)

title(main="A")

# second, go to IR_spectroscopy_functions.R and change dataset in fresh (line 19)
filt_ref_table <- read_reference_table("ref_table_filt_6_ID_sample_bottom.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("spectra_matrix_filt_6_ID_sample_bottom.txt", dir="el")

# ID_speciment, herbarium, top
best_PCs_value = 20
dapc2 <- dapc(filt_spect_matrix, filt_ref_table$ploidy, n.pca = best_PCs_value, n.da=1)

scatter(dapc2,
        cex = 0.5, 
        legend = TRUE, posi.leg = "topright", cleg = 0.8, txt.leg=c("2x","4x"),
        clabel = F,
        scree.pca = T, posi.pca = "topleft", 
        xax = 1, yax = 1, inset.solid = .5)

title(main="B")

# loading plot: which component differentiate more between groups

contribA <- loadingplot(dapc1$var.contr,axis = 1, thres = 1, 
                       lab.jitter = 10,
                       xlab = "wavelength (nm)",
                       #srt = 60,
                       main="loading plot")

contribB <- loadingplot(dapc2$var.contr,
                        #ylim=c(0,0.06),
                        axis = 1, thres = 1, 
                        lab.jitter = 10, 
                        xlab = "wavelength (nm)",
                        #srt = 60,
                        main="loading plot"
                        )

par(mfrow=c(1,1))



# Neural Network analysis ####
filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_individual_filt_6.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_individual_filt_6.txt", dir="el")

filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_sample_filt_6.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_sample_filt_6.txt", dir="el")

filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_individual_filt_6_bottom.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_individual_filt_6_bottom.txt", dir="el")

filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_individual_filt_6_top.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_individual_filt_6_top.txt", dir="el")

filt_ref_table <- read_reference_table("ref_table_filt_6.txt", dir="filt")
filt_spect_matrix <- read_spect_matrix("spectra_matrix_filt_6.txt", dir="filt")

results_nn_train_test <- train_test_nn(spect_matrix = filt_spect_matrix, 
                                                   ref_table = filt_ref_table,
                                                   seed = 222,
                                                   hidden_neurons = c(37, 25))

results_nn_LOO <- leave_one_out_nn(spect_matrix = filt_spect_matrix, 
                                   ref_table = filt_ref_table,
                                   seed = 222,
                                   hidden_neurons = c(36, 11, 4))

results_nn_k_fold <- k_fold_test_nn(spect_matrix = filt_spect_matrix, 
                                    ref_table = filt_ref_table,
                                    seed = 222,
                                    k_fold = 25,
                                    hidden_neurons = c(36, 11, 4),
                                    quiet = F
                                    )

# Herbarium results(accuracy) with 35, 18 hidden neurons ID_sample dataset:
# k=37 -> 0.78
# k=74 -> 0.79

# Herbarium results (accuracy) with (37, 25) hidden neurons ID_individual
# k = 5 -> 0.81
# k = 10 -> 0.70
# k = 20 -> 0.75
# k = 30 -> 0.74
# k = 35 -> 0.77
# k = 74 -> 0.81 # this is different from LOO test (0.70!)

# 
# from https://www.heatonresearch.com/2017/06/01/hidden-layers.html
# The number of hidden neurons should be between the size of the input layer and the size of the output layer.
# The number of hidden neurons should be 2/3 the size of the input layer, plus the size of the output layer.
# The number of hidden neurons should be less than twice the size of the input layer.

res_optimization <- find_optimal_n_neurons(spect_matrix = filt_spect_matrix, 
                       ref_table = filt_ref_table,
                       ratio_neuron_layers = 0.4,
                       layers_range = c(2:2),
                       neurons_range = c(20:90), 
                       step_neurons = 5,
                       method = "k-fold",
                       k_fold = 20,
                       )
res_optimization$accuracy


res_optimization <- find_optimal_n_neurons_multicore(cpus=3,
                                           spect_matrix = filt_spect_matrix, 
                                           ref_table = filt_ref_table,
                                           ratio_neuron_layers = 0.4,
                                           layers_range = c(2:2),
                                           neurons_range = c(20:60), 
                                           step_neurons = 5,
                                           method = "k-fold",
                                           k_fold = 20,
)

plot(10:80,res_optimization[[4]], col="red",type = "l", 
     xlab = "neurons on first hidden layer",
     ylab = "accuracy of prediction")
# max accuracy 0.8918919 with c(37, 25) neurons (~0.66 ration between neuron layers)

sink("best_neuron_structure.txt")
print(res_optimization)
sink()

optimal_dimensions <- find_optimal_nnet_dimensions(spect_matrix = filt_spect_matrix, 
                                                   ref_table = filt_ref_table,
                                                   method = "k-fold",
                                                   k_fold =20,
                                                   seed = 222,
                                                   step_neur_range = 0.1,
                                                   ratio_neur_range = c(0.4,0.7),
                                                   layers_range = c(2:2),
                                                   neurons_range = c(20:80), 
                                                   step_neurons = 5
                                                   )

optimal_dimensions <- find_optimal_nnet_dimensions_multicore(
                                                   cpus = 3,
                                                   spect_matrix = filt_spect_matrix, 
                                                   ref_table = filt_ref_table,
                                                   method = "k-fold",
                                                   k_fold =5,
                                                   seed = 222,
                                                   step_neur_range = 0.2,
                                                   ratio_neur_range = c(0.4,0.6),
                                                   layers_range = c(2:2),
                                                   neurons_range = c(40:60), 
                                                   step_neurons = 10
)
# with FRESH material (ID_ind):
# k=11
# 46,19: 54.54% accuracy
# 48, 21: 60.60%
# 80,36: 63.63%
# 124, 53: 63.63
# k=33
# 56, 23=54.54%
# 58,25=57.57% -> 39% with k=33
# 58,33=60.60% -> with k=33, 60.0 % accuracy
# 87,53 = 64.61% -> this reduces to 53% using k=33 (LOO)
# using ID_sample
# 




# with HERBARIUM material:
# best results with k=10:
# 38, 16 = 82.43% accuracy
# 32, 15 = 85.13% accuracy
# with k=37 (1/2 of dataset size)
# 34, 14 = 82.43% accuracy
# 39, 17 = 83.78%
# 35, 18 = 85.13 (ratio=0.49)



# best ratio of hidden layers found = 0.41: 57, 24 with 94,59% accuracy

# Leave one out
succ_LOO <- leave_one_out_nn(spect_matrix = filt_spect_matrix, 
                             ref_table=filt_ref_table, 
                             hidden_neurons = c(36,11,4),
                             seed = 222,
                             quiet=T)
succ_LOO
succ_LOO_2<-k_fold_test_nn(spect_matrix = filt_spect_matrix, 
                           ref_table=filt_ref_table, 
                           k_fold=nlevels(filt_ref_table$ID_ind),
                           hidden_neurons = c(36,11,4),
                           seed = 222,
                           quiet=T)
{ 
  beepr::beep()
  Sys.sleep(0.8)
  beepr::beep()
}



require(plot3D)

plot(accuracy_values)
accuracy_values <- res_optimization[[3]]
colnames(accuracy_values)
slice3D(accuracy_values,
       x=as.integer(rownames(accuracy_values)),
       y=as.integer(colnames(accuracy_values)),
       colvar = )
require(smoothAPC)
plot3d(accuracy_values)


# perform DAPC and NNET k-fold test on Veronica spicata or V. longifolia alone #####
# 6 DAPC and NNET in total: 1. V. spicata fresh and 2. herbarium, 3. V. longifolia fresh and 4. herbarium, 
# 5. V.spicata combined-fresh-dry, 6. V.longifolia combined-fresh-dry

# before running next line, in IR_spectroscopy_functions.R change spectra_dataset in "fresh" (line 19)
filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_individual_filt_6.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_individual_filt_6.txt", dir="el")

spi_spect_matrix <- filt_spect_matrix[filt_ref_table$species=="spicata",]
spi_ref_table <- filt_ref_table[filt_ref_table$species=="spicata",]
spi_ref_table$ID_ind <- as.factor(as.character(spi_ref_table$ID_ind))
spi_ref_table$ID_spect<-as.factor(as.character(spi_ref_table$ID_spect))

lon_spect_matrix <- filt_spect_matrix[filt_ref_table$species=="longifolia",]
lon_ref_table <- filt_ref_table[filt_ref_table$species=="longifolia",]
lon_ref_table$ID_ind <- as.factor(as.character(lon_ref_table$ID_ind))
lon_ref_table$ID_spect<-as.factor(as.character(lon_ref_table$ID_spect))

best_PC_spi = DAPC_find_best_PCs(spi_ref_table, spi_spect_matrix, 10, 120, 10, length(unique(spi_ref_table$ID_ind)))
dapc_spi <- dapc(spi_spect_matrix, spi_ref_table$ploidy, n.pca = best_PC_spi, n.da=1)

best_PC_lon=DAPC_find_best_PCs(lon_ref_table, lon_spect_matrix,10,120,10,length(unique(lon_ref_table$ID_ind)))
dapc_lon <- dapc(lon_spect_matrix, lon_ref_table$ploidy, n.pca = best_PC_lon, n.da=1)

par(mfrow=c(2,2)) # to have 2 plots in one image
scatter(dapc_spi,
        cex = 0.5, 
        legend = TRUE, posi.leg = "topright", cleg = 0.8, txt.leg=c("2x","4x"),
        clabel = F,
        scree.pca = T, posi.pca = "topleft", 
        xax = 1, yax = 1, inset.solid = .5)

scatter(dapc_lon,
        cex = 0.5, 
        legend = TRUE, posi.leg = "topright", cleg = 0.8, txt.leg=c("2x","4x"),
        clabel = F,
        scree.pca = T, posi.pca = "topleft", 
        xax = 1, yax = 1, inset.solid = .5)

contribA <- loadingplot(dapc_spi$var.contr,axis = 1, thres = 1, 
                        xlab = "wavelength (nm)",
                        main="loading plot")

contribB <- loadingplot(dapc_lon$var.contr,
                        axis = 1, thres = 1, 
                        xlab = "wavelength (nm)",
                        main="loading plot"
)

par(mfrow=c(1,1)) # to have 2 plots in one image

# check results DAPC Vs NNET
best_PC = 20

results <- k_fold_test_nn_and_DAPC(hidden_neurons = c(36,11,4),
                                   n_PCA = best_PC,
                                   spect_matrix=lon_spect_matrix, 
                                   ref_table=lon_ref_table, 
                                   k_fold = nrow(lon_spect_matrix), #use leave-one-out since we don't have much samples in this case
                                   seed = 222,
                                   quiet = F)

succ_NN=sum(results[[1]]$pred_success_NN==1)
succ_DAPC=sum(results[[1]]$pred_success_DAPC==1)
n=nlevels(lon_ref_table$ID_ind)

prop.test(x=c(succ_NN,succ_DAPC),
          n=c(n,n))

# RESULTS:
# 2-sample test for equality of proportions with continuity correction
# 
# data:  c(succ_NN, succ_DAPC) out of c(n, n)
# X-squared = 3.75, df = 1, p-value = 0.05281
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   -0.77217153 -0.02782847
# sample estimates:
#   prop 1    prop 2 
# 0.4666667 0.8666667 

results <- k_fold_test_nn_and_DAPC(hidden_neurons = c(36,11,4),
                                   n_PCA = best_PC,
                                   spect_matrix=spi_spect_matrix, 
                                   ref_table=spi_ref_table, 
                                   k_fold = nrow(spi_spect_matrix),
                                   seed = 222,
                                   quiet = F)


succ_NN=sum(results[[1]]$pred_success_NN==1)
succ_DAPC=sum(results[[1]]$pred_success_DAPC==1)
n=nlevels(spi_ref_table$ID_ind)

prop.test(x=c(succ_NN,succ_DAPC),
          n=c(n,n))

# RESULTS
# 2-sample test for equality of proportions with continuity correction
# 
# data:  c(succ_NN, succ_DAPC) out of c(n, n)
# X-squared = 0.45714, df = 1, p-value = 0.499
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   -0.5396793  0.2063460
# sample estimates:
#   prop 1    prop 2 
# 0.3333333 0.5000000 

# before running next line, in IR_spectroscopy_functions.R change spectra_dataset in "dry" (line 19)
filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_individual_filt_6.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_individual_filt_6.txt", dir="el")

spi_spect_matrix <- filt_spect_matrix[filt_ref_table$species=="spicata",]
spi_ref_table <- filt_ref_table[filt_ref_table$species=="spicata",]
spi_ref_table$ID_ind <- as.factor(as.character(spi_ref_table$ID_ind))
spi_ref_table$ID_spect<-as.factor(as.character(spi_ref_table$ID_spect))

lon_spect_matrix <- filt_spect_matrix[filt_ref_table$species=="longifolia",]
lon_ref_table <- filt_ref_table[filt_ref_table$species=="longifolia",]
lon_ref_table$ID_ind <- as.factor(as.character(lon_ref_table$ID_ind))
lon_ref_table$ID_spect<-as.factor(as.character(lon_ref_table$ID_spect))

best_PC_spi=DAPC_find_best_PCs(spi_ref_table,spi_spect_matrix, 10, 120, 10, length(unique(spi_ref_table$ID_ind)))
dapc3 <- dapc(spi_spect_matrix, spi_ref_table$ploidy, n.pca = best_PC_spi, n.da=1)

best_PC_lon=DAPC_find_best_PCs(lon_ref_table,lon_spect_matrix, 10, 120, 10, length(unique(lon_ref_table$ID_ind)))
dapc4 <- dapc(lon_spect_matrix, lon_ref_table$ploidy, n.pca = best_PC_lon, n.da=1)

# plot all results in a unique graph
par(mfrow=c(2,2),
    mai=c(.7,.7,.2,.1)#A numerical vector of the form c(bottom, left, top, right) which gives the margin size specified in inches.
)
#par(mfrow=c(2,2))
{
  scatter(dapc1,
          cex = 0.5, 
          legend = TRUE, posi.leg = "top", cleg = 0.8, txt.leg=c("2x","4x"),
          clabel = F,
          #scree.pca = T, posi.pca = "topleft", 
          xax = 1, yax = 1, inset.solid = .5)
  #title(main="A")
  scatter(dapc2,
          cex = 0.5, 
          legend = TRUE, posi.leg = "top", cleg = 0.8, txt.leg=c("2x","4x"),
          clabel = F,
          #scree.pca = T, posi.pca = "topleft", 
          xax = 1, yax = 1, inset.solid = .5)
  #title(main="B")
  
  scatter(dapc3,
          cex = 0.5, 
          legend = TRUE, posi.leg = "top", cleg = 0.8, txt.leg=c("2x","4x"),
          clabel = F,
          #scree.pca = T, posi.pca = "topleft", 
          xax = 1, yax = 1, inset.solid = .5)
  #title(main="C")
  scatter(dapc4,
          cex = 0.5, 
          legend = TRUE, posi.leg = "top", cleg = 0.8, txt.leg=c("2x","4x"),
          clabel = F,
          #scree.pca = T, posi.pca = "topleft", 
          xax = 1, yax = 1, inset.solid = .5)
  #title(main="D")
}
par(mfrow=c(1,1),mai=c(1,1,1,1))

# check results DAPC Vs NNET
best_PC = 20

results <- k_fold_test_nn_and_DAPC(hidden_neurons = c(36,11,4),
                                   n_PCA = best_PC,
                                   spect_matrix=lon_spect_matrix, 
                                   ref_table=lon_ref_table, 
                                   k_fold = nrow(lon_spect_matrix), #use leave-one-out since we don't have much samples in this case
                                   seed = 222,
                                   quiet = F)

succ_NN=sum(results[[1]]$pred_success_NN==1)
succ_DAPC=sum(results[[1]]$pred_success_DAPC==1)
n=nlevels(lon_ref_table$ID_ind)

prop.test(x=c(succ_NN,succ_DAPC),
          n=c(n,n))

# RESULTS
# 2-sample test for equality of proportions with continuity correction
# 
# data:  c(succ_NN, succ_DAPC) out of c(n, n)
# X-squared = 0.58217, df = 1, p-value = 0.4455
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   -0.1257354  0.3419516
# sample estimates:
#   prop 1    prop 2 
# 0.7567568 0.6486486 


results <- k_fold_test_nn_and_DAPC(hidden_neurons = c(36,11,4),
                                   n_PCA = best_PC,
                                   spect_matrix=spi_spect_matrix, 
                                   ref_table=spi_ref_table, 
                                   k_fold = nrow(spi_spect_matrix),
                                   seed = 222,
                                   quiet = F)


succ_NN=sum(results[[1]]$pred_success_NN==1)
succ_DAPC=sum(results[[1]]$pred_success_DAPC==1)
n=nlevels(lon_ref_table$ID_ind)

prop.test(x=c(succ_NN,succ_DAPC),
          n=c(n,n))
# RESULTS
# 2-sample test for equality of proportions with continuity correction
# 
# data:  c(succ_NN, succ_DAPC) out of c(n, n)
# X-squared = 1.133, df = 1, p-value = 0.2871
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   -0.08855848  0.35882875
# sample estimates:
#   prop 1    prop 2 
# 0.8108108 0.6756757 


#combined dataset
# before running next line, in IR_spectroscopy_functions.R change spectra_dataset in "combined_dry_fresh" (line 21)
filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_individual_filt_6.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_individual_filt_6.txt", dir="el")

spi_spect_matrix <- filt_spect_matrix[filt_ref_table$species=="spicata",]
spi_ref_table <- filt_ref_table[filt_ref_table$species=="spicata",]
spi_ref_table$ID_ind <- as.factor(as.character(spi_ref_table$ID_ind))
spi_ref_table$ID_spect<-as.factor(as.character(spi_ref_table$ID_spect))

lon_spect_matrix <- filt_spect_matrix[filt_ref_table$species=="longifolia",]
lon_ref_table <- filt_ref_table[filt_ref_table$species=="longifolia",]
lon_ref_table$ID_ind <- as.factor(as.character(lon_ref_table$ID_ind))
lon_ref_table$ID_spect<-as.factor(as.character(lon_ref_table$ID_spect))

best_PC_spi=DAPC_find_best_PCs(spi_ref_table,spi_spect_matrix, 10, 120, 10, length(unique(spi_ref_table$ID_ind)))
best_PC_lon=DAPC_find_best_PCs(lon_ref_table,lon_spect_matrix, 10, 120, 10, length(unique(lon_ref_table$ID_ind)))

best_PC_lon = 40
best_PC_spi = 10

results <- k_fold_test_nn_and_DAPC(hidden_neurons = c(41,17,7),
                                   n_PCA = best_PC_lon,
                                   spect_matrix=lon_spect_matrix, 
                                   ref_table=lon_ref_table, 
                                   k_fold = nrow(lon_spect_matrix), #use leave-one-out since we don't have much samples in this case
                                   seed = 222,
                                   quiet = F)

succ_NN=sum(results[[1]]$pred_success_NN==1)
succ_DAPC=sum(results[[1]]$pred_success_DAPC==1)
n=nlevels(lon_ref_table$ID_ind)

prop.test(x=c(succ_NN,succ_DAPC),
          n=c(n,n))

# RESULTS
# 2-sample test for equality of proportions with continuity correction
# 
# data:  c(succ_NN, succ_DAPC) out of c(n, n)
# X-squared = 1.1712, df = 1, p-value = 0.2792
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   -0.30734006  0.07657083
# sample estimates:
#   prop 1    prop 2 
# 0.6538462 0.7692308 

results <- k_fold_test_nn_and_DAPC(hidden_neurons = c(41,17,7),
                                   n_PCA = best_PC_spi,
                                   spect_matrix=spi_spect_matrix, 
                                   ref_table=spi_ref_table, 
                                   k_fold = nrow(spi_spect_matrix),
                                   seed = 222,
                                   quiet = F)


succ_NN=sum(results[[1]]$pred_success_NN==1)
succ_DAPC=sum(results[[1]]$pred_success_DAPC==1)
n=nlevels(lon_ref_table$ID_ind)

prop.test(x=c(succ_NN,succ_DAPC),
          n=c(n,n))

# RESULTS

# 2-sample test for equality of proportions with continuity correction
# 
# data:  c(succ_NN, succ_DAPC) out of c(n, n)
# X-squared = 6.0352e-31, df = 1, p-value = 1
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   -0.2107903  0.1723288
# sample estimates:
#   prop 1    prop 2 
# 0.7115385 0.7307692 


# statistics on results ####

# ROSE package to correct for umbalanced dataset
library(ROSE)

# compute the average of the difference and test if statistically is different from 0
# perform t-test (corrected for correlation - t-test assumes observations to be independent, but here we have the overlapping dataset)
# see https://stats.stackexchange.com/questions/343830/how-to-test-to-test-if-one-model-is-better-than-another-on-independent-dataset
# paper at link https://jmlr.org/papers/v18/16-305.html
# repeat k fold m times for both DAPC and NNET, so to get a vector of differences of accuracies



# best results NNet: herbarium, ID individual, both leaf sides (neurons = 36, 11, 4, PC=60, k_fold=25)
# IMPORTANT!: change line 17 in IR_spectroscopy_functions.R and reload
filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_individual_filt_6.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_individual_filt_6.txt", dir="el")

results_1 <- k_fold_test_nn_and_DAPC(hidden_neurons = c(36,11,4),
                          n_PCA = 60,
                          spect_matrix=filt_spect_matrix, 
                          ref_table=filt_ref_table, 
                          k_fold = 25,
                          seed = 222,
                          quiet = F)
  
{ 
  beepr::beep()
  Sys.sleep(0.8)
  beepr::beep()
}
results_1<-results_1[[1]]
write.csv(results_1, file="results_herbarium_NN_36_11_4_Vs_DAPC_60PC_k25.csv")
setwd("~/Dropbox/Infrared Spectroscopy/R script analysis/NN_Vs_DAPC_results/")

results_1 <- read.csv("results_herbarium_NN_36_11_4_Vs_DAPC_60PC_k25.csv", header=T)
results_1$real_ploidy<-as.factor(results_1$real_ploidy)
results_1$pred_ploidy_NN <- as.factor(results_1$pred_ploidy_NN)
results_1$pred_ploidy_DAPC <- as.factor(results_1$pred_ploidy_DAPC)


# McNemar's test
#contingency table
contingengy_table_1 <- table(results_1$pred_success_NN,results_1$pred_success_DAPC)
mcnemar.test(contingengy_table_1, correct = T)
# McNemar's Chi-squared test with continuity correction
# 
# data:  contingengy_table_1
# McNemar's chi-squared = 1.8947, df = 1, p-value = 0.1687


# Z-test
succ_NN=sum(results_1$pred_success_NN==1)
succ_DAPC=sum(results_1$pred_success_DAPC==1)
n=nlevels(filt_ref_table$ID_ind)

prop.test(x=c(succ_NN,succ_DAPC),
          n=c(n,n))

# RESULTS
# 2-sample test for equality of proportions with continuity correction
# 
# data:  c(succ_NN, succ_DAPC) out of c(74, 74)
# X-squared = 1.5439, df = 1, p-value = 0.214
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   -0.04590399  0.23509318
# sample estimates:
#   prop 1    prop 2 
# 0.8513514 0.7567568 

# generate new balanced data by ROSE

table(results_1$real_ploidy)
balanced_results_1 <- ROSE(real_ploidy ~ ., data=results_1, seed=123)$data
balanced_results_1$pred_success_DAPC = ifelse(balanced_results_1$real_ploidy==balanced_results_1$pred_ploidy_DAPC,1,0)
balanced_results_1$pred_success_NN = ifelse(balanced_results_1$real_ploidy==balanced_results_1$pred_ploidy_NN,1,0)
table(balanced_results_1$real_ploidy)

View(results_1)
View(balanced_results_1)
# Z-test on balanced data-set

succ_NN=sum(balanced_results_1$pred_success_NN==1)
succ_DAPC=sum(balanced_results_1$pred_success_DAPC==1)
n=nrow(balanced_results_1)

prop.test(x=c(succ_NN,succ_DAPC),
          n=c(n,n))
# RESULTS
# 2-sample test for equality of proportions with continuity correction
# 
# data:  c(succ_NN, succ_DAPC) out of c(n, n)
# X-squared = 0.14968, df = 1, p-value = 0.6988
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   -0.1097345  0.1908156
# sample estimates:
#   prop 1    prop 2 
# 0.7837838 0.7432432 

# best results NNet: combined, ID individual, both leaf sides (neurons = 41, 17, 7, 50PC, k_fold=35)
# IMPORTANT!: change line 17 in IR_spectroscopy_functions.R and reload
filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_individual_filt_6.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_individual_filt_6.txt", dir="el")

results_2 <- k_fold_test_nn_and_DAPC(hidden_neurons = c(41,17,7),
                                 n_PCA = 50,
                                 spect_matrix=filt_spect_matrix, 
                                 ref_table=filt_ref_table, 
                                 k_fold = 35,
                                 seed = 222,
                                 quiet = F)

{ 
  beepr::beep()
  Sys.sleep(0.8)
  beepr::beep()
}
results_2<-results_2[[1]]
write.csv(results_2, file="results_combined_NN_41_17_7_Vs_DAPC_50PC_k25.csv")

results_2 <- read.csv("results_combined_NN_41_17_7_Vs_DAPC_50PC_k25.csv", header=T)
results_2$real_ploidy <- as.factor(results_2$real_ploidy)
results_2$pred_ploidy_NN <- as.factor(results_2$pred_ploidy_NN)
results_2$pred_ploidy_DAPC <- as.factor(results_2$pred_ploidy_DAPC)
# McNemar's test
#contingency table
contingengy_table_2 <- table(results_2$pred_success_NN, results_2$pred_success_DAPC)
mcnemar.test(contingengy_table_2, correct = T)
# McNemar's Chi-squared test with continuity correction
# 
# data:  contingengy_table_2
# McNemar's chi-squared = 2.3824, df = 1, p-value = 0.1227

# Z-test
succ_NN=sum(results_2$pred_success_NN==1)
succ_DAPC=sum(results_2$pred_success_DAPC==1)
n=nlevels(filt_ref_table$ID_ind)
prop.test(x=c(succ_NN,succ_DAPC),
          n=c(n,n))
# 2-sample test for equality of proportions with continuity
# correction
# 
# data:  c(succ_NN, succ_DAPC) out of c(n, n)
# X-squared = 1.8056, df = 1, p-value = 0.179
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   -0.03793264  0.22484852
# sample estimates:
#   prop 1    prop 2 
# 0.7476636 0.6542056 


# generate new balanced data by ROSE

table(results_2$real_ploidy)
balanced_results_2 <- ROSE(real_ploidy ~ ., data=results_2, seed=123)$data
balanced_results_2$pred_success_DAPC = ifelse(balanced_results_2$real_ploidy==balanced_results_2$pred_ploidy_DAPC,1,0)
balanced_results_2$pred_success_NN = ifelse(balanced_results_2$real_ploidy==balanced_results_2$pred_ploidy_NN,1,0)
table(balanced_results_2$real_ploidy)

View(results_2)
View(balanced_results_2)
# Z-test on balanced data-set

succ_NN=sum(balanced_results_2$pred_success_NN==1)
succ_DAPC=sum(balanced_results_2$pred_success_DAPC==1)
n=nrow(balanced_results_2)

prop.test(x=c(succ_NN,succ_DAPC),
          n=c(n,n))
# RESULTS
# 
# 2-sample test for equality of proportions with continuity correction
# 
# data:  c(succ_NN, succ_DAPC) out of c(n, n)
# X-squared = 0.88096, df = 1, p-value = 0.3479
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   -0.06069145  0.19153257
# sample estimates:
#   prop 1    prop 2 
# 0.7757009 0.7102804 

# best results DAPC: herbarium, ID individual, top leaf side (PC = 50) Vs NNet (24,8,3), k_fold=25
# IMPORTANT!: change line 17 in IR_spectroscopy_functions.R and reload
filt_ref_table <- read_reference_table("ref_table_average_spectra_ID_individual_filt_6_top.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("average_spectra_ID_individual_filt_6_top.txt", dir="el")

results_3 <- k_fold_test_nn_and_DAPC(hidden_neurons = c(24,8,3),
                                 n_PCA = 50,
                                 spect_matrix=filt_spect_matrix, 
                                 ref_table=filt_ref_table, 
                                 k_fold = 25,
                                 seed = 222,
                                 quiet = F)

{ 
  beepr::beep()
  Sys.sleep(0.8)
  beepr::beep()
}
results_3 <- results_3[[1]]
write.csv(results_3, file="results_herbarium_NN_24_8_3_Vs_DAPC_50PC_k25.csv")

results_3 <- read.csv("results_herbarium_NN_24_8_3_Vs_DAPC_50PC_k25.csv", header=T)
results_3$real_ploidy <- as.factor(results_3$real_ploidy)
results_3$pred_ploidy_NN <- as.factor(results_3$pred_ploidy_NN)
results_3$pred_ploidy_DAPC <- as.factor(results_3$pred_ploidy_DAPC)

# McNemar's test
#contingency table
contingengy_table_3 <- table(results_3$pred_success_NN, results_3$pred_success_DAPC)
mcnemar.test(contingengy_table_3, correct = T)
# 
# McNemar's Chi-squared test with continuity correction
# 
# data:  contingengy_table_3
# McNemar's chi-squared = 2.0833, df = 1, p-value = 0.1489

# Z-test
succ_NN=sum(results_3$pred_success_NN==1)
succ_DAPC=sum(results_3$pred_success_DAPC==1)
n=nlevels(filt_ref_table$ID_ind)
prop.test(x=c(succ_NN,succ_DAPC),
          n=c(n,n))

#RESULTS
# 2-sample test for equality of proportions with continuity correction
# 
# data:  c(succ_NN, succ_DAPC) out of c(n, n)
# X-squared = 1.261, df = 1, p-value = 0.2615
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   -0.22709316  0.05318011
# sample estimates:
#   prop 1    prop 2 
# 0.7826087 0.8695652 

# generate new balanced data by ROSE

table(results_3$real_ploidy)
balanced_results_3 <- ROSE(real_ploidy ~ ., data=results_3, seed=123)$data
balanced_results_3$pred_success_DAPC = ifelse(balanced_results_3$real_ploidy==balanced_results_3$pred_ploidy_DAPC,1,0)
balanced_results_3$pred_success_NN = ifelse(balanced_results_3$real_ploidy==balanced_results_3$pred_ploidy_NN,1,0)
table(balanced_results_3$real_ploidy)

View(results_3)
View(balanced_results_3)
# Z-test on balanced data-set

succ_NN = sum(balanced_results_3$pred_success_NN == 1)
succ_DAPC = sum(balanced_results_3$pred_success_DAPC == 1)
n = nrow(balanced_results_3)

prop.test(x = c(succ_NN, succ_DAPC),
          n = c(n,n))
# RESULTS

# 2-sample test for equality of proportions with continuity correction
# 
# data:  c(succ_NN, succ_DAPC) out of c(n, n)
# X-squared = 0.97656, df = 1, p-value = 0.323
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   -0.20129578  0.05636824
# sample estimates:
#   prop 1    prop 2 
# 0.8260870 0.8985507 

# best results DAPC: fresh, ID sample, bottom leaf side (PCs = 20, 96,68,48 neurons, k_fold=11)
# IMPORTANT!: change line 17 in IR_spectroscopy_functions.R and reload
filt_ref_table <- read_reference_table("ref_table_filt_6_ID_sample_bottom.txt", dir="el")
filt_spect_matrix <- read_spect_matrix("spectra_matrix_filt_6_ID_sample_bottom.txt", dir="el")

results_4 <- k_fold_test_nn_and_DAPC(hidden_neurons = c(96,68,48),
                                 n_PCA = 20,
                                 spect_matrix=filt_spect_matrix, 
                                 ref_table=filt_ref_table, 
                                 k_fold = 11,
                                 seed = 222,
                                 quiet = F)

for (i in 1:5){ 
  beepr::beep()
  Sys.sleep(0.2)
}
results_4 <- results_4[[1]]
write.csv(results_4, file="results_fresh_NN_96_68_48_Vs_DAPC_20PC_k11.csv")

results_4 <- read.csv("results_fresh_NN_96_68_48_Vs_DAPC_20PC_k11.csv", header=T)
results_4$real_ploidy <- as.factor(results_4$real_ploidy)
results_4$pred_ploidy_NN <- as.factor(results_4$pred_ploidy_NN)
results_4$pred_ploidy_DAPC <- as.factor(results_4$pred_ploidy_DAPC)


# McNemar's test
#contingency table
contingengy_table_4 <- table(results_4[[1]]$pred_success_NN, results_4[[1]]$pred_success_DAPC)
mcnemar.test(contingengy_table_4, correct = T)
# 
# McNemar's Chi-squared test with continuity correction
# 
# data:  contingengy_table_4
# McNemar's chi-squared = 0.12121, df = 1, p-value = 0.7277

# Z-test

succ_NN=sum(results[[1]]$pred_success_NN==1)
succ_DAPC=sum(results[[1]]$pred_success_DAPC==1)
n=nrow(filt_ref_table)
prop.test(x=c(succ_NN,succ_DAPC),
          n=c(n,n))
# 2-sample test for equality of proportions with continuity correction
# 
# data:  c(succ_NN, succ_DAPC) out of c(n, n)
# X-squared = 0.13595, df = 1, p-value = 0.7123
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   -0.2249066  0.1325989
# sample estimates:
#   prop 1    prop 2 
# 0.6307692 0.6769231 

# generate new balanced data by ROSE

table(results_4$real_ploidy)
balanced_results_4 <- ROSE(real_ploidy ~ ., data = results_4, seed = 123)$data
balanced_results_4$pred_success_DAPC = ifelse(balanced_results_4$real_ploidy == balanced_results_4$pred_ploidy_DAPC, 1, 0)
balanced_results_4$pred_success_NN = ifelse(balanced_results_4$real_ploidy == balanced_results_4$pred_ploidy_NN, 1, 0)
table(balanced_results_4$real_ploidy)

View(results_4)
View(balanced_results_4)
# Z-test on balanced data-set

succ_NN = sum(balanced_results_4$pred_success_NN == 1)
succ_DAPC = sum(balanced_results_4$pred_success_DAPC == 1)
n = nrow(balanced_results_4)

prop.test(x = c(succ_NN, succ_DAPC),
          n = c(n,n))
# RESULTS
# 
# 2-sample test for equality of proportions with continuity correction
# 
# data:  c(succ_NN, succ_DAPC) out of c(n, n)
# X-squared = 1.1791, df = 1, p-value = 0.2775
# alternative hypothesis: two.sided
# 95 percent confidence interval:
#   -0.28865624  0.07327163
# sample estimates:
#   prop 1    prop 2 
# 0.5692308 0.6769231 
