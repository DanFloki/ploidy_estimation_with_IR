# FUNCTIONS utilized for IR spectral analysis
# AUTHOR: Daniele Buono
# daniele.buono@tu-dresden.de
# daniele.buo@gmail.com
# last edit: 31 May 2022

# required libraries
suppressMessages(library(lubridate)) # for execution time calculation
suppressMessages(library(neuralnet)) # for neural network
suppressMessages(library(signal)) # for S-Golay filter
suppressMessages(library(adegenet)) # for DAPC
suppressMessages(library(parallel)) # to exploid multicore processors and HPC



# constant values
wvl <- 350:2500

# NOTE: use folder "dry" for herbarium material, "fresh" for living plants, "comb" for combined dataset
 spectra_dataset <-"dry"
# spectra_dataset <- "combined_dry_fresh"
# spectra_dataset <- "fresh"

filt_data_dir <- paste0("/home/daniele/Dropbox/Infrared Spectroscopy/R script analysis/",spectra_dataset,"/filtered_data/")
elaborated_data_dir <- paste0("~/Dropbox/Infrared Spectroscopy/R script analysis/",spectra_dataset,"/elaborated_data/")
analysis_graphs_dir <- paste0("~/Dropbox/Infrared Spectroscopy/R script analysis/",spectra_dataset,"/analysis_graphs/")
analysis_result_tables_dir <- paste0("~/Dropbox/Infrared Spectroscopy/R script analysis/",spectra_dataset,"/result_tables/")


# read the spectra from several output files (.sed) and create a single matrix to store them
read_raw_spectrophotometer_data <- function(reference_table){
  #import spectra (sed-files)
  #load sed-file and locate position of spectrum values following the header data
  spectra_names <- paste0(reference_table$ID_spect,".sed")
  temp <- readLines(spectra_names[1])
  firstline <- which(substring(temp,1,3)=="Wvl") #exclude first lines (header with info) searching for Wvl
  test <- read.table(spectra_names[1],skip=firstline-1,header=TRUE,sep="\t")
  
  #generate spectrum matrix from all spectra included in reference value table
  spectra_matrix <- matrix(nrow=length(spectra_names),ncol=dim(test)[1])
  colnames(spectra_matrix) <- test[,1]
  rownames(spectra_matrix) <- rownames(reference_table)
  for (i in 1:length(spectra_names)){
    spectrum <- read.table(spectra_names[i],skip=firstline-1,header=TRUE,sep="\t")
    spectra_matrix[i,] <- spectrum[,2]
  }
  
  row.names(spectra_matrix) <- reference_table$ID_spect

  return(spectra_matrix)
}

# read matrix of spectra saved previously
read_spect_matrix <- function(file_name, dir){
  if(dir=="filt"){
    spect_matrix <- as.matrix(read.table(paste0(filt_data_dir, file_name), header=T, sep="\t"))
  } else if (dir=="el"){
    spect_matrix <- as.matrix(read.table(paste0(elaborated_data_dir, file_name), header=T, sep="\t"))
  } else {
    cat("ERROR: dir option =",dir,", but can be \"filt\" for filtered data or \"el\" for elaborated data!")
  }
  #remove letter X in front of each column name
  colnames(spect_matrix) <- substring(colnames(spect_matrix),2)
  return(spect_matrix)
}

# read reference table for matrix fo spectra
read_reference_table <- function (file_name, dir){
  if(dir=="filt"){
    ref_table <- read.table(paste0(filt_data_dir, file_name), header=T, sep="\t")
  } else if (dir=="el"){
    ref_table <- read.table(paste0(elaborated_data_dir, file_name), header=T, sep="\t")
  } else {
    cat("ERROR: dir option can be \"filt\" for filtered data or \"el\" for elaborated data!\n")
  }
  return(ref_table)
}

# remove spectra that shows > 100 % reflectance
remove_anomalous_spectra <- function(spectra_matrix){
  anomalous_spectra <- as.character()
  for (i in row.names(spectra_matrix)){ #for each spectra saved 
    for (j in colnames(spectra_matrix)){ # for each wavelength
      if(spectra_matrix[i,j]>100){
        cat ("spectra",i,"has reflectance > 100% at",j,"nm (",spectra_matrix[i,as.character(j)],")\n")
        anomalous_spectra<-c(anomalous_spectra, i)
        break
      }
    }
  }
  
  filt_spect_matrix <- spectra_matrix[!rownames(spectra_matrix) %in% anomalous_spectra, ]
  return(filt_spect_matrix)
}

# function to remove samples represented by less than 3 spectra
remove_sampled_less_than_3 <- function(ref_table, spect_matrix){
  filt_spect_matrix<-spect_matrix
  filt_ref_table<-ref_table
  for(i in levels(ref_table$ID_sample)){ # for each ID_sample
    if(sum(ref_table$ID_sample==i)<3){ # if there are less than 3 spectra
      cat(paste0("removing ",i,", n spectra = ",sum(ref_table$ID_sample==i),"\n"))
      filt_spect_matrix<-filt_spect_matrix[!rownames(filt_spect_matrix)%in%ref_table[ref_table$ID_sample==i,]$ID_spect,]
    }
  }
  return(filt_spect_matrix)
}

# function to apply Savitzky-Golay filter

apply_sgolay_filter <- function (spect_matrix, p, n){
  # p=3 and n=33 were chosen after inspecting some spectra, as below
  # # install.packages("signal")
  # library("signal")
  # 
  # #plot one spectra to have an idea on the power of the filter
  # spc2<-sgolayfilt(spect.matrix[200,], #x
  #                  3,# p=3, filter order
  #                  33,# n=p+3-p%%2,
  #                  0# m=0,
  #                  # ts=1
  # )
  # plot(spect.matrix[200,],type="l")
  # plot(spc2,type="l")
  
  filt_spect_matrix <- spect_matrix
  for (i in 1:nrow(filt_spect_matrix)){ #make it row by row because method can give back only a vector not a matrix
    filt_spect_matrix[i,]<-sgolayfilt(spect_matrix[i,],
                                             p, # filter order
                                             n) # n
  }
  return(filt_spect_matrix)
}

# function to remove spectra too similar to background
remove_similar_to_background <- function(spect_matrix, ref_table, spect_matrix_background, k){
  # make a plot to follow progress
  plot(wvl,spect_matrix_background, type="l", col="red",
       ylab="Reflect. %", xlab="wavelength", ylim=c(0,100), 
       main=paste0("Removed spectra, k=",k))
  lines(wvl,spect_matrix_background+k, type="l", lty=2, col="red")
  lines(wvl,spect_matrix_background-k, type="l", lty=2, col="red")
  
  index_remove<-numeric()
  for (i in 1:nrow(spect_matrix)){ 
    
    diff_i <- mean(abs(spect_matrix[i,]-spect_matrix_background))
    
    if(diff_i < k){
      cat(paste0("removing spectra: ", i,", ", as.character(ref_table[i,]$ID_spect),"\n"))
      lines(wvl, spect_matrix[i,], type="l", col="blue")
      index_remove<-c(index_remove, i)
    }
  }
  cat("removed",length(index_remove),"spectra.\n")
  
  lines(wvl,spect_matrix_background, type="l", lty=2, col="red")
  lines(wvl,spect_matrix_background+k, type="l", lty=2, col="red")
  lines(wvl,spect_matrix_background-k, type="l", lty=2, col="red")
  
  filt_spect_matrix <- spect_matrix[-index_remove,]
  return(filt_spect_matrix)
}

# reduce dataset size by averaging reflectance values on ranges of X wavelength
reduce_dataset_size <- function(spect_matrix, reduction_factor){
  # spect_matrix<-filt_5_spect_matrix
  # reduction_factor=5
  filt_spect_matrix<-matrix(nrow = nrow(spect_matrix), ncol = ncol(spect_matrix)/reduction_factor)
  dim(filt_spect_matrix)
  wvl_names<-as.character()
  for (i in seq(from=min(wvl), to=max(wvl)-1, by=reduction_factor)){
    wvl_names<-c(wvl_names, paste0(as.character(i),"-",as.character(i+reduction_factor-1)))
  }
  wvl_names
  # length(wvl_names)
  
  colnames(filt_spect_matrix) <- wvl_names
  row.names(filt_spect_matrix) <- row.names(spect_matrix)
  # e <- edit(filt_spect_matrix)
  
  for(i in 1:nrow(filt_spect_matrix)){
    for(j in 1:ncol(filt_spect_matrix)){
      k=(j*reduction_factor)-(reduction_factor-1)
      filt_spect_matrix[i,j]<-mean(spect_matrix[i, k:(k+(reduction_factor-1))])
    }
  }
  return(filt_spect_matrix)
}

# standard error calculation
std_err <- function(x) sd(x)/sqrt(length(x))


# plot spectra for each individual with std_err
plot_spectra_std_err <- function (ref_table, spect_matrix, col_graphs, max_y, avg_std_err_plot=T){
  par(mfrow=c(1,1))
  by_val = length(wvl)/ ncol(spect_matrix) #introduced to overcome plotting of reduced datasets
  wvl_number <- seq(min(wvl), max(wvl), by=by_val) #number of columns in reduced dataset (avg every 5 wvl)
  titles=c("A","B","C","D")
  for (i in unique(ref_table$ID_ind)){ #for each ID
    spectra_i = spect_matrix[i == ref_table$ID_ind,]
    ref_t_i = ref_table[ref_table$ID_ind == i, ]
    
    file_name_i = paste0(as.character(i),".svg")
    svg(file_name_i, width = 3.5, height = 3)

    
    if( !is.null(nrow(spectra_i)) ){
      plot(wvl_number, spectra_i[1, ], type="l", col = col_graphs,
           ylab = "Reflectance %", xlab="wavelength (nm)", ylim = c(0, max_y), 
           main =  titles[which(unique(ref_table$ID_ind)==i)])
 
      for (j in 2 : nrow(spectra_i)){
        lines(wvl_number, spectra_i[j, ], type="l", col = col_graphs)
      }
      if(avg_std_err_plot){
        #avg
        lines(wvl_number, colMeans(spectra_i), type="l", col="blue")
        #std error
        stdErr <- apply(spectra_i, 2, std_err)
        lines(wvl_number, colMeans(spectra_i) + stdErr, type="l", lty=2, col="blue")
        lines(wvl_number, colMeans(spectra_i) - stdErr, type="l", lty=2, col="blue")
      }
      
    } else {
      plot(wvl_number, spectra_t_o, type="l", col=col_graphs,
           ylab="Reflect. %", xlab="wavelength (nm)", ylim=c(0,max_y), 
           main=paste0("Spectra sample ",as.character(i),"_t_o, ",ref_t_i$species[1]))
    }
    dev.off()
    
  }
}

# plot avg_spectra only for each individual with std_err
plot_avg_spectra_std_err <- function (ref_table, spect_matrix, col_graphs, max_y){
  par(mfrow=c(1,1))
  
  by_val = length(wvl)/ ncol(spect_matrix) #introduced to overcome plotting of reduced datasets
  wvl_number <- seq(min(wvl), max(wvl), by=by_val) #number of columns in reduced dataset (avg every 5 wvl)
  titles=c("A","B","C","D")
  for (i in unique(ref_table$ID_ind)){ #for each ID
    spectra_i = spect_matrix[i == ref_table$ID_ind,]
    file_name_i = paste0("AVG_",as.character(i),".svg")
    svg(file_name_i, width = 3.5, height = 3)
    
    
    if( !is.null(nrow(spectra_i)) ){
      #avg
      plot(wvl_number, colMeans(spectra_i), type="l", col = col_graphs,
           ylab = "Reflectance %", xlab="wavelength (nm)", ylim = c(0, max_y), 
           main = titles[which(unique(ref_table$ID_ind)==i)])
      #std error
      stdErr <- apply(spectra_i, 2, std_err)
      lines(wvl_number, colMeans(spectra_i) + stdErr, type="l", lty=2, col=col_graphs)
      lines(wvl_number, colMeans(spectra_i) - stdErr, type="l", lty=2, col=col_graphs)
    } else {
      plot(wvl_number, spectra_t_o, type="l", col=col_graphs,
           ylab="Reflect. %", xlab="wavelength (nm)", ylim=c(0,max_y), 
           main=titles[which(unique(ref_table$ID_ind)==i)])
    }
    dev.off()
    
  }
}

# plot avg_spectra only for each individual with std_err
plot_avg_spectra_std_err_one_graph <- function (ref_table, spect_matrix, col_graphs = c("#1B9E77","#D95F02"), max_y, file_name) {
  par(mfrow=c(1,1))
  
  by_val = length(wvl)/ ncol(spect_matrix) #introduced to overcome plotting of reduced datasets
  wvl_number <- seq(min(wvl), max(wvl), by=by_val) #number of columns in reduced dataset (avg every 5 wvl)
  svg(paste0(file_name,".svg"), width = 10, height = 6)
  first=T
  
  for (i in unique(ref_table$ID_ind)){ #for each ID
    spectra_i = spect_matrix[i == ref_table$ID_ind,]

    if( !is.null(nrow(spectra_i)) ){
      #avg
      if(first){
        plot(wvl_number, colMeans(spectra_i), type="l", col = col_graphs,
             ylab = "Reflectance %", xlab="wavelength (nm)", ylim = c(0, max_y))
        first=F
      } else {
          lines(wvl_number, colMeans(spectra_i), type="l", col = col_graphs[which(unique(ref_table$ID_ind)==i)])
      }
      #std error
      stdErr <- apply(spectra_i, 2, std_err)
      lines(wvl_number, colMeans(spectra_i) + stdErr, type="l", lty=2, col=col_graphs[which(unique(ref_table$ID_ind)==i)])
      lines(wvl_number, colMeans(spectra_i) - stdErr, type="l", lty=2, col=col_graphs[which(unique(ref_table$ID_ind)==i)])
    } else {
      plot(wvl_number, spectra_t_o, type="l", col=col_graphs,
           ylab="Reflect. %", xlab="wavelength (nm)", ylim=c(0,max_y), 
           main=titles[which(unique(ref_table$ID_ind)==i)])
    }
  }
  dev.off()
}

# plot spectra in different files according to sample and sampling position
plot_spectra <- function (ref_table, spect_matrix, col_graphs, max_y){
  par(mfrow=c(1,1))

  by_val=length(wvl)/ncol(spect_matrix) #introduced to overcome plotting of reduced datasets
  wvl_number<-seq(min(wvl),max(wvl),by=by_val) #number of columns in reduced dataset (avg every 5 wvl)
  # length(wvl_number)
  # max(wvl_number)
  for (i in levels(ref_table$ID_ind)){ #for each ID
    
    spectra_i=spect_matrix[i == ref_table$ID_ind,]
    ref_t_i=ref_table[ref_table$ID_ind==i,]
    
    # T_O
    if(T %in% (ref_t_i$leaf_side=="t" & ref_t_i$leaf_age=="o")){
      file_name_i=paste0(as.character(i),"_t_o.svg")
      svg(file_name_i, width = 10, height = 6)
      spectra_t_o<-spectra_i[ref_t_i$leaf_side=="t" & ref_t_i$leaf_age=="o",]
      if( !is.null(nrow(spectra_t_o)) ){
        plot(wvl_number, spectra_t_o[1,], type="l", col=col_graphs,
             ylab="Reflectance %", xlab="wavelength", ylim=c(0,max_y), 
             main=paste0("Spectra sample ",as.character(i),"_t_o, ",ref_t_i$species[1]))
        for (j in 2:nrow(spectra_t_o)){
          lines(wvl_number, spectra_t_o[j,],type="l", col=col_graphs)
        }
        #avg
        lines(wvl_number, colMeans(spectra_t_o),type="l", col="blue")
        #std error
        stdErr<-apply(spectra_t_o, 2, std_err)
        lines(wvl_number,colMeans(spectra_t_o)+stdErr,type="l", lty=2, col="blue")
        lines(wvl_number,colMeans(spectra_t_o)-stdErr,type="l", lty=2, col="blue")
      } else {
        plot(wvl_number, spectra_t_o, type="l", col=col_graphs,
             ylab="Reflect. %", xlab="wavelength", ylim=c(0,max_y), 
             main=paste0("Spectra sample ",as.character(i),"_t_o, ",ref_t_i$species[1]))
      }
      dev.off()
    }
    
    #T_N
    if(T %in% (ref_t_i$leaf_side=="t" & ref_t_i$leaf_age=="n")){
      file_name_i=paste0(as.character(i),"_t_n.svg")
      svg(file_name_i, width = 10, height = 6)
      spectra_t_n<-spectra_i[ref_t_i$leaf_side=="t" & ref_t_i$leaf_age=="n",]
      if( !is.null(nrow(spectra_t_n)) ){
        plot(wvl_number, spectra_t_n[1,], type="l", col=col_graphs,
             ylab="Reflectance %", xlab="wavelength", ylim=c(0,max_y), 
             main=paste0("Spectra sample ",as.character(i),"_t_n, ",ref_t_i$species[1]))
        for (j in 2:nrow(spectra_t_n)){
          lines(wvl_number, spectra_t_n[j,],type="l", col=col_graphs)
        }
        #avg
        lines(wvl_number, colMeans(spectra_t_n),type="l", col="blue")
        #std error
        stdErr<-apply(spectra_t_n, 2, std_err)
        lines(wvl_number,colMeans(spectra_t_n)+stdErr,type="l", lty=2, col="blue")
        lines(wvl_number,colMeans(spectra_t_n)-stdErr,type="l", lty=2, col="blue")
      } else {
        plot(wvl_number, spectra_t_n, type="l", col=col_graphs,
             ylab="Reflectance %", xlab="wavelength", ylim=c(0,max_y), 
             main=paste0("Spectra sample ",as.character(i),"_t_n, ",ref_t_i$species[1]))
      }
      dev.off()
    }
    
    #B_O
    if(T %in% (ref_t_i$leaf_side=="b" & ref_t_i$leaf_age=="o")){
      file_name_i=paste0(as.character(i),"_b_o.svg")
      svg(file_name_i, width = 10, height = 6)
      spectra_b_o<-spectra_i[ref_t_i$leaf_side=="b" & ref_t_i$leaf_age=="o",]
      if( !is.null(nrow(spectra_b_o)) ){
        plot(wvl_number, spectra_b_o[1,], type="l", col=col_graphs,
             ylab="Reflectance %", xlab="wavelength", ylim=c(0,max_y), 
             main=paste0("Spectra sample ",as.character(i),"_b_o, ",ref_t_i$species[1]))
        for (j in 2:nrow(spectra_b_o)){
          lines(wvl_number, spectra_b_o[j,],type="l", col=col_graphs)
        }
        #avg
        lines(wvl_number, colMeans(spectra_b_o),type="l", col="blue")
        #std error
        stdErr<-apply(spectra_b_o, 2, std_err)
        lines(wvl_number,colMeans(spectra_b_o)+stdErr,type="l", lty=2, col="blue")
        lines(wvl_number,colMeans(spectra_b_o)-stdErr,type="l", lty=2, col="blue")
      } else {
        plot(wvl_number, spectra_b_o, type="l", col=col_graphs,
             ylab="Reflect. %", xlab="wavelength", ylim=c(0,max_y), 
             main=paste0("Spectra sample ",as.character(i),"_b_o, ",ref_t_i$species[1]))
      }
      dev.off()
    }
    
    #B_N
    if(T %in% (ref_t_i$leaf_side=="b" & ref_t_i$leaf_age=="n")){
      file_name_i=paste0(as.character(i),"_b_n.svg")
      svg(file_name_i, width = 10, height = 6)
      spectra_b_n<-spectra_i[ref_t_i$leaf_side=="b" & ref_t_i$leaf_age=="n",]
      if( !is.null(nrow(spectra_b_n)) ){
        plot(wvl_number, spectra_b_n[1,], type="l", col=col_graphs,
             ylab="Reflect. %", xlab="wavelength", ylim=c(0,max_y), 
             main=paste0("Spectra sample ",as.character(i),"_b_n, ",ref_t_i$species[1]))
        for (j in 2:nrow(spectra_b_n)){
          lines(wvl_number, spectra_b_n[j,],type="l", col=col_graphs)
        }
        #avg
        lines(wvl_number, colMeans(spectra_b_n),type="l", col="blue")
        #std error
        stdErr<-apply(spectra_b_n, 2, std_err)
        lines(wvl_number,colMeans(spectra_b_n)+stdErr,type="l", lty=2, col="blue")
        lines(wvl_number,colMeans(spectra_b_n)-stdErr,type="l", lty=2, col="blue")
      } else {
        plot(wvl_number, spectra_b_n, type="l", col=col_graphs,
             ylab="Reflect. %", xlab="wavelength", ylim=c(0,max_y), 
             main=paste0("Spectra sample ",as.character(i),"_b_n, ",ref_t_i$species[1]))
      }
      dev.off()
    }
  }
}

plot_ind_spectra <- function (spect_matrix, col_graphs, max_y){
    
    by_val=length(wvl)/ncol(spect_matrix) #introduced to overcome plotting of reduced datasets
    wvl_number<-seq(min(wvl),max(wvl),by=by_val) #number of columns in reduced dataset (avg every 5 wvl)

    for (i in row.names(spect_matrix)){ #for each ID
        file_name_i=paste0(as.character(i),".svg")
        svg(file_name_i, width = 10, height = 6)

        plot(wvl_number, spect_matrix[i,], type="l", col=col_graphs,
               ylab="Reflect. %", xlab="Wavelength (nm)", ylim=c(0,max_y), 
               main=paste0("Spectra ", as.character(i)))
        dev.off()
    }
}

# average spectra of samples
average_sample_spectra <- function (ref_table, spect_matrix){
  sample_spect_matrix <- matrix(nrow=nlevels(ref_table$ID_sample),
                              ncol=ncol(spect_matrix)) # a column for every wavelength
  row.names(sample_spect_matrix) = levels(ref_table$ID_sample) # give a name to the rows (2n and 4n)
  colnames(sample_spect_matrix) = colnames(spect_matrix) #give a name to the columns (from 350 to 2500)
  dim(sample_spect_matrix)
  
  ref_avg <- data.frame(row.names = c("ID_spect", "species", "ploidy", "leaf_age", "leaf_side",	"date",	"age_side", "ID_ind",	"ID_sample"))
  for(i in row.names(sample_spect_matrix)){ # for every ID_sample
    for(j in colnames(sample_spect_matrix)){ #for every wavelenght
      sample_spect_matrix[i,j] <- mean(spect_matrix[ref_table$ID_sample == i, j])
    }
    ref_avg <- rbind(ref_avg, head(ref_table[ref_table$ID_sample == i, -1], 1)) # create a new reference table adding only one line of the information on 6 spectra
  }
  ref_avg$ID_spect=ref_avg$ID_sample
  sample_spect_ref <- list("ref_table" = ref_avg, "spect_matrix" = sample_spect_matrix)
  return(sample_spect_ref)
}


# average spectra of individuals
average_individual_spectra <- function (ref_table, 
                                        spect_matrix){
  
  individual_spect_matrix <- matrix(nrow = nlevels(ref_table$ID_ind),
                                ncol = ncol(spect_matrix)) # a column for every wavelength
  row.names(individual_spect_matrix) = levels(ref_table$ID_ind) # give a name to the rows (2n and 4n)
  colnames(individual_spect_matrix) = colnames(spect_matrix) #give a name to the columns (from 350 to 2500)
  #dim(individual_spect_matrix)
  
  ref_avg <-data.frame()
   for(i in row.names(individual_spect_matrix)){ # for every ID_sample
    for(j in colnames(individual_spect_matrix)){ #for every wavelenght
      individual_spect_matrix[i,j] <- mean(spect_matrix[ref_table$ID_ind==i,j])
    }
    ref_avg <- rbind(ref_avg, head(ref_table[ref_table$ID_ind==i,-1],1)) #create a new reference table adding only one line of the informations on 6 spectra
  }
  
  # reducing information on reference table (ex. ID_spectra)
  ref_avg <- ref_avg[, c(1,2,3,6,8)]
  ref_avg$ID_spect <- ref_avg$ID_ind
  sample_spect_ref <- list("ref_table" = ref_avg, "spect_matrix" = individual_spect_matrix)
  return(sample_spect_ref)
}

DAPC_find_best_PCs <- function(ref_table, spect_matrix, min_PCA, max_PCA, step_increase, rep_i){
  

  # 
  # # # just for debug____________________________________________________
  # ref_table = filt_ref_table
  # spect_matrix = filt_spect_matrix
  # min_PCA = 10
  # max_PCA = 200
  # step_increase = 10
  # rep_i = nlevels(ref_table$ID_ind) # to test with all individuals
  # rep_i=10
  # # __________________________________________________________________
  # 
  
  perc_succ_nPCA = numeric()
  time0 = Sys.time()
  best_PC_all = 0
  highest_pred_success = 0
  for (n_PCA in seq(min_PCA, max_PCA, by = step_increase)){
    
    # ___ for final statistics __
    pred_success <- numeric()
    ID_spectra <- character()
    species <- character()
    ploidy <- character()
    pred_ploidy <- character()
    post_P <- numeric()
    #____________________________
    
    sampled = character() # to avoid sampling the same individual more than once
    
    for(i in 1:rep_i){ # repeat rep_i times:
      cat("\rcalculating error using ", n_PCA, " PCs, spectra ", i,"/", rep_i, "... ")
      # extract randomly one individual from dataset, avoiding repeats
      repeat{
        out_ID <- sample(levels(ref_table$ID_ind), 1)
        if(!(as.character(out_ID) %in% sampled)){break}
      }
      sampled = c(sampled, as.character(out_ID))
      
      excluded_ID_spectra <- ref_table[ref_table$ID_ind == out_ID, ]$ID_spect
      j <- which(rownames(spect_matrix) %in% excluded_ID_spectra) # get index of spectra to exclude
      train_spect <- spect_matrix[-j, ] # exclude spectra
      train_sp <- ref_table$ploidy[-j] # exclude info on ploidy level of excluded sample
      
      # construct DAPC model with rest of samples
      train_dapc <- dapc(train_spect, train_sp, n.pca=n_PCA, n.da=1)
      
      # predict ploidy for every spectra leaved out and record success
      excluded_spectra <- spect_matrix[rownames(spect_matrix) %in% as.character(excluded_ID_spectra), ]
      ploidy_current_i <- ref_table[ref_table$ID_spect %in% excluded_ID_spectra, ]$ploidy
      if( !is.null(length(ploidy_current_i)) ){ ploidy_current_i=ploidy_current_i[1]} #in case there are more lines indicating ploidy level
      
      # collect info for statistics ____________________________________________
      ID_spectra <- c(ID_spectra, as.character(excluded_ID_spectra))
      species <- c(species, as.character(ref_table[ref_table$ID_ind == out_ID, ]$species))
      ploidy <- c(ploidy, as.character(ref_table[ref_table$ID_ind == out_ID, ]$ploidy))
      
      #_________________________________________________________________________
      if( !is.null(nrow(excluded_spectra)) ){ # nrow gives back NULL when there is only one line in matrix (since it is not a matrix anymore but a vector)
        for (k in 1:nrow(excluded_spectra)){ 
          test_spect <- as.data.frame(t(excluded_spectra[k,]))
          pred_dapc <- predict.dapc(train_dapc,newdata=test_spect)
          
          pred_success <- c(pred_success, ifelse(as.character(pred_dapc$assign) == ploidy_current_i,1,0))
          pred_ploidy <- c(pred_ploidy, as.character(pred_dapc$assign))
          post_P <- c(post_P, ifelse(pred_dapc$assign == "2", pred_dapc$posterior[1], pred_dapc$posterior[2]))
        }
      } else {
        test_spect <- as.data.frame(t(excluded_spectra))
        pred_dapc <- predict.dapc(train_dapc, newdata = test_spect)

        pred_success <- c(pred_success, ifelse(as.character(pred_dapc$assign) == ploidy_current_i,1,0))
        pred_ploidy <- c(pred_ploidy, as.character(pred_dapc$assign))
        post_P <- c(post_P, ifelse(pred_dapc$assign == "2", pred_dapc$posterior[1], pred_dapc$posterior[2]))
        
      }
    }
    
    # _collect statistics in one file each number of PCs tested _____
    test_stats <- data.frame(ID_spectra, species, ploidy, pred_ploidy, pred_success, post_P)
    
    file_name <- paste0(analysis_result_tables_dir,"stats_LOO_",n_PCA,"PCs.csv")
    write.table(test_stats, file_name, row.names=F, sep="\t")
    # ________________________________________________________________
    
    perc_succ_current_PCA = sum(pred_success == 1) / length(pred_success)
    perc_succ_nPCA = c(perc_succ_nPCA, perc_succ_current_PCA)
    
    if(perc_succ_current_PCA > highest_pred_success){
      best_PC_all = n_PCA
      highest_pred_success = perc_succ_current_PCA
      cat("best PC value till now =",best_PC_all,"\n")
    }
    
    #time left estimation
    time_i = round(as.numeric(as.period(Sys.time()-time0, unit="sec")),0)
    perc_compl <- n_PCA/max_PCA
    t_tot = time_i/perc_compl
    s_left = t_tot-time_i
    m_left = as.integer(s_left / 60)
    h_left = as.integer(m_left / 60)
    s_left = as.integer(s_left %% 60)
    m_left = as.integer(m_left %% 60)
    cat("DONE. Accuracy = ", sum(pred_success==1) / length(pred_success),".\n")
    cat("Advanced ",paste0(round(perc_compl*100,2),"%")," time left approx.:",paste0(h_left,"h:",m_left,"m:",s_left,"s"),"\n")
    if(perc_compl == 1){
      beepr::beep()
      Sys.sleep(0.8)
      beepr::beep()
    }
  } #loop for each PCA
  
  
  
  #best_PC_all <- (which.max(perc_succ_nPCA)*step_increase) + min_PCA
  highest_score <- max(perc_succ_nPCA)
  
  plot(seq(min_PCA,max_PCA,by=step_increase),perc_succ_nPCA, ylim=c(0,1), type="b", lty=2,
       ylab="prediction success", xlab="PCs retained", 
       main=paste0("DAPC prediction success changing number of PCs\nY = ploidy rep_i = ", rep_i,
                   ", best PCs = ", best_PC_all,", score = ", round(highest_score, 3)))
  lines(seq(min_PCA,max_PCA,by=step_increase), rep(highest_score, length(seq(min_PCA,max_PCA,by=step_increase))), col="red")
  
  svg_file_name <- paste0(analysis_graphs_dir,"results xval LATEST ",Sys.time(),".svg")
  svg(file = svg_file_name, width = 8, height = 5)
  plot(seq(min_PCA,max_PCA,by=step_increase),perc_succ_nPCA, ylim=c(0,1), type="b", lty=2,
       ylab="prediction success", xlab="PCs retained", 
       main=paste0("DAPC prediction success changing number of PCs\nY = ploidy rep_i = ", rep_i,
                   ", best PCs = ", best_PC_all,", score = ", round(highest_score, 3)))
  lines(seq(min_PCA,max_PCA,by=step_increase), rep(highest_score, length(seq(min_PCA,max_PCA,by=step_increase))), col="red")
  dev.off()
  return(best_PC_all)
}

# function to assess performances of neural network by splitting dataset in two: training and testing
train_test_nn <- function (hidden_neurons = c(8,6), 
                           ref_table, 
                           spect_matrix, 
                           seed = 222, 
                           quiet = F){
 
  #normalization input data from 0 to 1
  spect_matrix = spect_matrix/100
  
  # dividing dataset into training and test by sampling ID_individuals randomly
  test_samples <- sample(levels(ref_table$ID_ind), size = nlevels(ref_table$ID_ind)/2)
  train_samples <- ref_table[!(ref_table$ID_ind %in% test_samples),]$ID_ind
  
  test_data_frame <- data.frame("ID_spectra" = row.names(spect_matrix[ref_table$ID_ind %in% test_samples,]), 
                                "ploidy" = ref_table[ref_table$ID_ind %in% test_samples,]$ploidy,
                                "species" = ref_table[ref_table$ID_ind %in% test_samples,]$species,
                                spect_matrix[ref_table$ID_ind %in% test_samples,])

  train_data_frame <- data.frame("ID_spectra" = row.names(spect_matrix[ref_table$ID_ind %in% train_samples,]), 
                                "ploidy" = ref_table[ref_table$ID_ind %in% train_samples,]$ploidy, 
                                "species" = ref_table[ref_table$ID_ind %in% train_samples,]$species,
                                spect_matrix[ref_table$ID_ind %in% train_samples,])
  
  # converting ploidy variable in 0 for 2x and 1 for 4x
  train_data_frame$ploidy = ifelse(train_data_frame$ploidy == 2, 0, 1)
  # converting ploidy variable in 0 for 2x and 1 for 4x
  test_data_frame$ploidy = ifelse(test_data_frame$ploidy == 2, 0, 1)
  
  # create formula to use in neural network ("ploidy ~ X350.354 + X355.359 + ... + X2495.2499")
  formula = paste("ploidy ~", paste(names(train_data_frame)[-c(1, 2, 3)], collapse = " + "))
  #cat("neural networks will be build by using the formula:", formula,"\n\n")
  set.seed(seed)
  # build model
  nn = neuralnet(formula,
                 data = train_data_frame, 
                 hidden = hidden_neurons, #number of hidden neurons in each layer
                 algorithm = "rprop+", # resilient backpropagation with weight backtracking
                 #learningrate =0.0001,
                 act.fct = "logistic",
                 err.fct = "sse", #ce or sse (cross-entropy or sum of square errors)
                 linear.output = FALSE,
                 lifesign = ifelse(!quiet,"full","none"),
                 rep = 1,
                 stepmax = 1000000,
                 threshold = 0.01) # threshold for partial derivates of the error function
  # plot(nn) # this take a loooong time because of the many inputs (X350, X351, ...) to the first node of the network
  
  # model prediction accuracy
  pred_train <- compute(nn, train_data_frame[, -c(1,2,3)])
  # extracts results of predictions and converts in binary
  p_train <- pred_train$net.result
  pred_train <- ifelse(p_train > 0.5, 4, 2)
  
  # make predictions of ploidy based on the TEST dataset minus ploidy level
  pred_test <- compute(nn, test_data_frame[, -c(1,2,3)])
  #extracts results of predictions and converts in binary
  p_test <- pred_test$net.result
  pred_test <- ifelse(p_test > 0.5, 4, 2)
  
  # converting ploidy variable in 0 for 2x and 1 for 4x
  train_data_frame$ploidy = ifelse(train_data_frame$ploidy == 0, 2, 4)
  # converting ploidy variable in 0 for 2x and 1 for 4x
  test_data_frame$ploidy = ifelse(test_data_frame$ploidy == 0, 2, 4)
  
  # create confusion table
  tab_train <- table("predicted" = pred_train, "real" = train_data_frame$ploidy)
  accuracy_train <- sum(diag(tab_train))/sum(tab_train)

  
  # create confusion table
  tab_test <- table(pred_test, test_data_frame$ploidy)
  accuracy_test <- sum(diag(tab_test))/sum(tab_test)
  
  results = list("TEST_samples" = test_samples, "TRAIN_samples" = train_samples,
                 "predictions_TRAIN","real" = pred_train,"predictions_TEST" = pred_test,
                 "conf_table_TRAIN" = tab_train,"conf_table_TEST" = tab_test,
                 "accuracy" = accuracy_test, "accuracy_train" = accuracy_train,
                 "neurons" = hidden_neurons)
  
  if(!quiet){
    cat("prediction accuracy based on TRAIN dataset =", accuracy_train, "; confusion table:\n")
    print(tab_train)
    cat("prediction accuracy based on TEST dataset =",accuracy_test,"; confusion table:\n")
    print(tab_test)
    file_name <- paste0(analysis_result_tables_dir,"NeuralNet/","stats_train-test_NNET_",paste(as.character(hidden_neurons), collapse = "_"),"_neurons_", Sys.time(),".csv")
    sink(file_name)
    print(results)
    sink()
  }

  return(results)
}


# k_fold_test DAPC
k_fold_test_DAPC <- function(spect_matrix, 
                             ref_table, 
                             k_fold = 10,
                             seed = 222,
                             quiet = F,
                             n_PCA=20){

  #normalization input data from 0 to 1
  
  # ___ for final statistics __
  pred_success <- numeric()
  ID_spectra <- character()
  species <- character()
  ploidy <- character()
  pred_ploidy <- character()
  post_P <- numeric()
  #____________________________
  
  # shuffle dataset order
  set.seed(seed)
  shuffled_IDs = sample(as.character(levels(ref_table$ID_ind)))
  
  # split dataset into k groups
  f = numeric()
  k_i = 1
  for(i in 1:length(shuffled_IDs)){
    f = c(f, k_i)
    k_i = k_i+1
    if(k_i == k_fold+1) {k_i = 1}
  }
  groups = split(shuffled_IDs, f)
  if(!quiet){
    cat("starting k-fold cross validation test\nGroups:")
    print(groups)
  }
  
  # for each unique group
  for(k_i in 1 : k_fold){ #repeat k times:
    time0 = Sys.time()
    
    # training group (all samples that are not in the group k_i)
    excluded_ID_spectra <- ref_table[ref_table$ID_ind %in% groups[[k_i]],]$ID_spect
    j <- which(rownames(spect_matrix) %in% excluded_ID_spectra) # get index of spectra to exclude
    train_spect <- spect_matrix[-j, ] # exclude spectra
    train_sp <- ref_table$ploidy[-j] # exclude info on ploidy level of excluded sample
    
    # construct DAPC model with rest of samples
    train_dapc <- dapc(train_spect, train_sp, n.pca=n_PCA, n.da=1)
    
    # predict ploidy for every spectra leaved out and record success
    excluded_spectra <- spect_matrix[rownames(spect_matrix) %in% as.character(excluded_ID_spectra), ]
    ploidy_current_i <- ref_table[ref_table$ID_spect %in% excluded_ID_spectra, ]$ploidy
    if( !is.null(length(ploidy_current_i)) ){ ploidy_current_i=ploidy_current_i[1]} #in case there are more lines indicating ploidy level
    
    if( !is.null(nrow(excluded_spectra)) ){ # nrow gives back NULL when there is only one line in matrix (since it is not a matrix anymore but a vector)
      for (j in 1:nrow(excluded_spectra)){ 
        test_spect <- as.data.frame(t(excluded_spectra[j,]))
        pred_dapc <- predict.dapc(train_dapc, newdata=test_spect)
        
        pred_success <- c(pred_success, ifelse(as.character(pred_dapc$assign) == ploidy_current_i,1,0))
        pred_ploidy <- c(pred_ploidy, as.character(pred_dapc$assign))
        post_P <- c(post_P, ifelse(pred_dapc$assign == "2", pred_dapc$posterior[1], pred_dapc$posterior[2]))
      }
    } else {
      test_spect <- as.data.frame(t(excluded_spectra))
      pred_dapc <- predict.dapc(train_dapc, newdata = test_spect)
      
      pred_success <- c(pred_success, ifelse(as.character(pred_dapc$assign) == ploidy_current_i,1,0))
      pred_ploidy <- c(pred_ploidy, as.character(pred_dapc$assign))
      post_P <- c(post_P, ifelse(pred_dapc$assign == "2", pred_dapc$posterior[1], pred_dapc$posterior[2]))
      
    }
    # collect info for statistics ____________________________________________
    ID_spectra <- c(ID_spectra, as.character(excluded_ID_spectra))
    species <- c(species, as.character(ref_table[ref_table$ID_ind == out_ID, ]$species))
    ploidy <- c(ploidy, as.character(ref_table[ref_table$ID_ind == out_ID, ]$ploidy))

   
    #time left estimation
    if(!quiet){
      time_i = round(as.numeric(as.period(Sys.time()-time0, unit="sec")),0)
      perc_compl = (k_i/k_fold)*100
      t_tot = time_i*k_fold
      s_left = t_tot-(t_tot*(perc_compl/100))
      m_left = as.integer(s_left/60)
      h_left = as.integer(m_left/60)
      s_left = as.integer(s_left%%60)
      m_left = as.integer(m_left%%60)
      cat("\radvanced ",paste0(round(perc_compl,2),"%")," time left approx.:",
          paste0(h_left,"h:",m_left,"m:",s_left,"s"))
    }
    
  }
  
  test_stats <- data.frame(ID_spectra, species, ploidy, pred_ploidy, pred_success, post_P)
  
  accuracy = sum(pred_success == 1)/length(pred_success)
  
  if(!quiet){
    print("Done. Accuracy:",accuracy)
  }
  
  results = list(test_stats, "accuracy" = accuracy)
  return(results)
}

# k-fold test neural networks
k_fold_test_nn <- function(hidden_neurons = c(37,25),
                           spect_matrix, 
                           ref_table, 
                           k_fold = 10,
                           seed = 222,
                           quiet = F){
  
  set.seed(seed)
  #normalization input data from 0 to 1
  spect_matrix = spect_matrix/100
  
  # for final statistics_______
  pred_success <- numeric()
  real_ploidy <- character()
  pred_ploidy <- character()
  prob <- numeric()
  sampled <- character()
  #____________________________
  
  data_frame <- data.frame("ID_spectra" = ref_table$ID_spect, 
                           "ploidy" = ifelse(ref_table$ploidy == 2, 0, 1),   # converting ploidy variable in 0 for 2x and 1 for 4x
                           "species" = ref_table$species,
                           spect_matrix)
  data_frame$ID_spectra <- as.character(data_frame$ID_spectra)
  
  # create formula to use in neural network ("ploidy ~ X350.354 + X355.359 + ... + X2495.2499")
  formula <- paste("ploidy ~", paste(names(data_frame[,-c(1,2,3)]), collapse = " + "))

  # shuffle dataset order
  shuffled_IDs = sample(as.character(levels(ref_table$ID_ind)))

  # split dataset into k groups
  f = numeric()
  k_i = 1
  for(i in 1:length(shuffled_IDs)){
    f = c(f, k_i)
    k_i = k_i+1
    if(k_i == k_fold+1) {k_i = 1}
  }
  groups = split(shuffled_IDs, f)
  if(!quiet){
    cat("starting k-fold cross validation test\nGroups:")
    print(groups)
  }
  
  # for each unique group
  for(k_i in 1 : k_fold){ #repeat k times:
    time0 = Sys.time()

    #construct neuralnetwork with all samples that are not in the group k_i
    
    excluded_ID_spectra <- ref_table[ref_table$ID_ind %in% groups[[k_i]],]$ID_spect
    train_data_frame = data_frame[!rownames(data_frame) %in% excluded_ID_spectra, ]
    test_data_frame = data_frame[rownames(data_frame) %in% excluded_ID_spectra, ]
    
    train_nn = neuralnet(formula,
                   data = train_data_frame,
                   hidden = hidden_neurons, #number of hidden neurons in each layer
                   algorithm = "rprop+", # resilient backpropagation with weight backtracking
                   #learningrate =0.0001,
                   act.fct = "logistic",
                   err.fct = "sse", #ce or sse (cross-entropy or sum of square errors)
                   linear.output = FALSE,
                   lifesign = "none",
                   rep = 1,
                   stepmax = 1000000,
                   threshold = 0.01) # threshold for partial derivates of the error function
    
    if(is.null(train_nn$net.result)){
      cat("\n NO CONVERGENGE ON NEURAL NETWORK!\n")
    }
    # make predictions of ploidy based on the TEST dataset minus ploidy level
    pred_test <- predict(train_nn, test_data_frame[, -c(1,2,3)])
    #extracts results of predictions and converts in binary
    
    prob = c(prob, pred_test)
    pred_test <- ifelse(pred_test > 0.5, 4, 2)
    
    # converting ploidy variable in 0 for 2x and 1 for 4x
    test_data_frame$ploidy = ifelse(test_data_frame$ploidy == 0, 2, 4)
    
    sampled = c(sampled, row.names(test_data_frame))
    pred_ploidy = c(pred_ploidy, pred_test)
    real_ploidy = c(real_ploidy, test_data_frame$ploidy)
    pred_success = c(pred_success, ifelse(pred_test == test_data_frame$ploidy,1,0))
    
    #time left estimation
    if(!quiet){
      time_i = round(as.numeric(as.period(Sys.time()-time0, unit="sec")),0)
      perc_compl = (k_i/k_fold)*100
      t_tot = time_i*k_fold
      s_left = t_tot-(t_tot*(perc_compl/100))
      m_left = as.integer(s_left/60)
      h_left = as.integer(m_left/60)
      s_left = as.integer(s_left%%60)
      m_left = as.integer(m_left%%60)
      cat("\radvanced ",paste0(round(perc_compl,2),"%")," time left approx.:",
          paste0(h_left,"h:",m_left,"m:",s_left,"s"))
    }
    
  }
  
  accuracy = sum(pred_success == 1)/length(pred_success)
  test_stats <- data.frame("ID_spectra" = sampled, 
                           real_ploidy, 
                           pred_ploidy, 
                           pred_success,
                           prob)
  
  if(!quiet){
    file_name <- paste0(analysis_result_tables_dir,"NeuralNet/","stats_k_",k_fold,"_fold_NNET_", paste(as.character(hidden_neurons), collapse = "_"),"_neurons_", Sys.time(),".csv")
    write.table(test_stats, file_name, row.names = F, sep="\t")
  }
  
  results = list(test_stats, "accuracy" = accuracy, "neurons" = hidden_neurons)
  return(results)
}


# k-fold test neural networks and DAPC together
k_fold_test_nn_and_DAPC <- function(hidden_neurons = c(37,25),
                           n_PCA = 40,
                           spect_matrix, 
                           ref_table, 
                           k_fold = 10,
                           seed = 222,
                           quiet = F){
  
  
  #normalization input data from 0 to 1
  norm_spect_matrix = spect_matrix/100
  
  # for final statistics_______
  pred_success_NN <- numeric()
  pred_success_DAPC <- numeric()
  real_ploidy <- character()
  pred_ploidy_NN <- character()
  pred_ploidy_DAPC <- character()
  prob_NN <- numeric()
  sampled <- character()
  #____________________________
  
  data_frame <- data.frame("ID_spectra" = ref_table$ID_spect, 
                           "ploidy" = ifelse(ref_table$ploidy == 2, 0, 1),   # converting ploidy variable in 0 for 2x and 1 for 4x
                           "species" = ref_table$species,
                           norm_spect_matrix)
  data_frame$ID_spectra <- as.character(data_frame$ID_spectra)
  
  # create formula to use in neural network ("ploidy ~ X350.354 + X355.359 + ... + X2495.2499")
  formula <- paste("ploidy ~", paste(names(data_frame[,-c(1,2,3)]), collapse = " + "))
  
  # shuffle dataset order
  set.seed(seed)
  shuffled_IDs = sample(as.character(levels(ref_table$ID_ind)))
  
  # split dataset into k groups
  f = numeric()
  k_i = 1
  for(i in 1:length(shuffled_IDs)){
    f = c(f, k_i)
    k_i = k_i+1
    if(k_i == k_fold+1) {k_i = 1}
  }
  groups = split(shuffled_IDs, f)
  if(!quiet){
    cat("starting k-fold cross validation test\nGroups:")
    print(groups)
  }
  
  # for each unique group
  for(k_i in 1 : k_fold){ #repeat k times:
    time0 = Sys.time()
    
    #construct neuralnetwork with all samples that are not in the group k_i
    
    excluded_ID_spectra <- ref_table[ref_table$ID_ind %in% groups[[k_i]],]$ID_spect
    train_data_frame = data_frame[!rownames(data_frame) %in% excluded_ID_spectra, ]
    train_spectra = spect_matrix[!rownames(spect_matrix) %in% excluded_ID_spectra,]
    test_data_frame = data_frame[rownames(data_frame) %in% excluded_ID_spectra, ]
    test_spectra = spect_matrix[rownames(spect_matrix) %in% excluded_ID_spectra,]
    
    train_nn = neuralnet(formula,
                         data = train_data_frame,
                         hidden = hidden_neurons, #number of hidden neurons in each layer
                         algorithm = "rprop+", # resilient backpropagation with weight backtracking
                         #learningrate =0.0001,
                         act.fct = "logistic",
                         err.fct = "sse", #ce or sse (cross-entropy or sum of square errors)
                         linear.output = FALSE,
                         lifesign = "none",
                         rep = 1,
                         stepmax = 1000000,
                         threshold = 0.01) # threshold for partial derivates of the error function
    
    if(is.null(train_nn$net.result)){
      cat("\n NO CONVERGENGE ON NEURAL NETWORK!\n")
    }
    # DAPC predictions
    train_dapc <- dapc(train_spectra, train_data_frame$ploidy, n.pca=n_PCA, n.da=1)
    
    
    # NNET predictons (make predictions of ploidy based on the TEST dataset minus ploidy level)
    pred_test_NN <- predict(train_nn, test_data_frame[, -c(1,2,3)])
    
    if(k_fold == nlevels(ref_table$ID_ind)){
      pred_test_DAPC <- predict.dapc(train_dapc, newdata=t(test_spectra))
      
    } else {
      pred_test_DAPC <- predict.dapc(train_dapc, newdata=test_spectra)
    }
    
    #extracts results of predictions and converts in binary
    prob_NN = c(prob_NN, pred_test_NN)
    pred_test_NN <- ifelse(pred_test_NN > 0.5, 4, 2)
    pred_test_DAPC <- ifelse(pred_test_DAPC$assign==0,2,4)
    
    # converting ploidy variable in 0 for 2x and 1 for 4x
    test_data_frame$ploidy = ifelse(test_data_frame$ploidy == 0, 2, 4)
    
    sampled = c(sampled, row.names(test_data_frame))
    pred_ploidy_NN = c(pred_ploidy_NN, pred_test_NN)
    pred_ploidy_DAPC = c(pred_ploidy_DAPC, pred_test_DAPC)
    real_ploidy = c(real_ploidy, test_data_frame$ploidy)
    pred_success_NN = c(pred_success_NN, ifelse(pred_test_NN == test_data_frame$ploidy,1,0))
    pred_success_DAPC = c(pred_success_DAPC, ifelse(pred_test_DAPC==test_data_frame$ploidy,1,0))
    #time left estimation
    if(!quiet){
      time_i = round(as.numeric(as.period(Sys.time()-time0, unit="sec")),0)
      perc_compl = (k_i/k_fold)*100
      t_tot = time_i*k_fold
      s_left = t_tot-(t_tot*(perc_compl/100))
      m_left = as.integer(s_left/60)
      h_left = as.integer(m_left/60)
      s_left = as.integer(s_left%%60)
      m_left = as.integer(m_left%%60)
      cat("\radvanced ",paste0(round(perc_compl,2),"%")," time left approx.:",
          paste0(h_left,"h:",m_left,"m:",s_left,"s"))
    }
    
  }
  
  accuracy_NN = sum(pred_success_NN == 1)/length(pred_success_NN)
  accuracy_DAPC = sum(pred_success_DAPC == 1)/length(pred_success_DAPC)
  
  test_stats <- data.frame("ID_spectra" = sampled, 
                           real_ploidy, 
                           pred_ploidy_NN,
                           prob_NN,
                           pred_ploidy_DAPC,
                           pred_success_NN,
                           pred_success_DAPC)
  
  if(!quiet){
    file_name <- paste0(analysis_result_tables_dir,"NeuralNet/","stats_k_",k_fold,"_fold_NNET_", paste(as.character(hidden_neurons), collapse = "_"),"_neurons_", Sys.time(),".csv")
    write.table(test_stats, file_name, row.names = F, sep="\t")
  }
  
  results = list(test_stats, "accuracy_NN" = accuracy_NN, "neurons" = hidden_neurons,"accuracy_DAPC" = accuracy_DAPC)
  return(results)
}

# function to find optimal layers of neurons and number of neurons for each layer based on results 
# of neural_network_train_test() function
find_optimal_n_neurons <- function (ref_table, 
                                    spect_matrix, 
                                    method="t-t",
                                    seed = 222, 
                                    layers_range = c(1:5), 
                                    neurons_range = c(1:20),
                                    ratio_neuron_layers = 0.66,
                                    step_neurons = 2,
                                    quiet = F,
                                    #rep_i = 10,
                                    k_fold = 10){
  best_accuracy = 0
  best_neuron_value = 0
  best_result = NULL
  neuron_values = seq(min(neurons_range), max(neurons_range), by = step_neurons)
  accuracy_variability = matrix(ncol = length(layers_range),
                                nrow = length(neuron_values))
  colnames(accuracy_variability) = as.character(layers_range)
  rownames(accuracy_variability) = as.character(neuron_values)
  
  
  for(i in layers_range){
    if(!quiet) {cat("\ntesting with", i, "layers...\n")}
    for(j in neuron_values){
      hidden_neurons_value = rep(j, i)
      if(i>1){
        for(k in 2:length(hidden_neurons_value)){
          hidden_neurons_value[k] = as.integer((hidden_neurons_value[k-1]*ratio_neuron_layers))+1
        }
      }
      if(method == "t-t"){
        results_ij <- train_test_nn(spect_matrix = spect_matrix, 
                                                ref_table = ref_table,
                                                seed = seed,
                                                hidden_neurons = hidden_neurons_value,
                                                quiet = T)
      } else if (method == "LOO") {
        results_ij <- leave_one_out_nn(spect_matrix = spect_matrix, 
                                       ref_table = ref_table,
                                       seed = seed,
                                       hidden_neurons = hidden_neurons_value,
                                       quiet = T)
      } else if (method == "k-fold"){
        results_ij <- k_fold_test_nn(spect_matrix = spect_matrix, 
                                   ref_table = ref_table, 
                                   k_fold = k_fold,
                                   hidden_neurons = hidden_neurons_value,
                                   seed = seed,
                                   quiet = T)
      }
      else {
        cat("ERROR: method test can be only \"t-t\" for train-test, \"LOO\" for leave-one-out, or \"k-fold\" for k-fold cross (value passed:", method,")!\n")
      }
      
      accuracy_variability[as.character(j), as.character(i)] = results_ij$accuracy
      if(results_ij$accuracy > best_accuracy){
        best_accuracy = results_ij$accuracy
        best_neuron_value = hidden_neurons_value
        best_result = results_ij
        if(!quiet) {cat("\rbest value for layers", best_neuron_value, "with", best_accuracy*100, "% accuracy")}
      }
    }
  }
  
  final_results <- list(best_accuracy, best_result, best_neuron_value, accuracy_variability)
  
  return(final_results)
}

find_optimal_n_neurons_multicore <- function (cpus = 3,
                                    ref_table, 
                                    spect_matrix, 
                                    method="t-t",
                                    seed = 222, 
                                    layers_range = c(1:5), 
                                    neurons_range = c(20:80),
                                    ratio_neuron_layers = 0.66,
                                    step_neurons = 2,
                                    quiet = F,
                                    k_fold = 10){
  best_accuracy = 0
  best_neuron_value = 0
  best_result = NULL
  neuron_values = seq(min(neurons_range), max(neurons_range), by = step_neurons)
  accuracy_variability = matrix(ncol = length(layers_range),
                                nrow = length(neuron_values))
  colnames(accuracy_variability) = as.character(layers_range)
  rownames(accuracy_variability) = as.character(neuron_values)
  
  
  for(i in layers_range){ # repeat for each number of layers
    if(!quiet) {cat("\ntesting with", i, "layers...\n")}
    
    # calculate number of neurons in each layer according to ratio and number of layers
    hidden_neurons_values = list()
    for(j in neuron_values){
      hidden_neurons_value_i = rep(j, i)
      if(i>1){
        for(k in 2:length(hidden_neurons_value_i)){
          hidden_neurons_value_i[k] = as.integer((hidden_neurons_value_i[k-1]*ratio_neuron_layers))+1
        }
      }
      hidden_neurons_values[[length(hidden_neurons_values)+1]] <- hidden_neurons_value_i
    }
    
    
    if(method == "t-t"){
      results <- mclapply(hidden_neurons_values,
                          train_test_nn, 
                          ref_table=ref_table, 
                          spect_matrix = spect_matrix,
                          seed = seed,
                          quiet = T, 
                          mc.cores = cpus)
      
    } else if (method == "LOO") {
      results <- mclapply(hidden_neurons_values,
                        leave_one_out_nn,
                        spect_matrix = spect_matrix,
                        ref_table = ref_table,
                        quiet=T,
                        seed = seed,
                        mc.cores = cpus)
      
    } else if (method == "k-fold"){
      results <- mclapply(hidden_neurons_values,
                          k_fold_test_nn,
                          spect_matrix = spect_matrix,
                          ref_table = ref_table, 
                          k_fold = k_fold,
                          seed = seed,
                          quiet = T,
                          mc.cores=cpus)
    }
    else {
      cat("ERROR: method test can be only \"t-t\" for train-test, \"LOO\" for leave-one-out, or \"k-fold\" for k-fold cross (value passed:", method,")!\n")
    }
  }
  
  for(res_i in results){
    if(res_i$accuracy>best_accuracy){
      best_accuracy = res_i$accuracy
      best_result = res_i
      best_neuron_value = res_i$neurons
    }
  }
  final_results <- list("best_accuracy"=best_accuracy, "best_results"= best_result, "best_neuron_value"=best_neuron_value)
  
  return(final_results)
}


# function to find optimal dimensions of neural network, including ratio between hidden layer sizes
find_optimal_nnet_dimensions <- function (ref_table, # reference table for the spectra
                                          spect_matrix, # matrix containing spectra
                                          seed = 222, # random seed
                                          method ="t-t",
                                          layers_range = c(1:3), # min and max number of hidden layers to test
                                          neurons_range = c(20:60), # min and max number of neurons on the first hidden layer to test
                                          step_neurons = 2, # step of increase from min number of neurons to test
                                          ratio_neur_range = c(0.1,1.0), # min and max value of ratio between hidden layers to test
                                          step_neur_range = 0.1,
                                          k_fold = 10,
                                          quiet=F) {
  if(!quiet){
    cat("looking for optiman nnet dimensions, starting: ", as.character(Sys.time()))
    cat("\nseed =",seed,", method =",method," k =",k_fold,",layers range =", layers_range,",neurons range = (", min(neurons_range),",",max(neurons_range),")\n")
    cat("step neurons =",step_neurons,",ratio neur range =", ratio_neur_range,",step neur range =", step_neur_range)
  }

  
  best_accuracy = 0
  best_ratio_layers = 0
  time0 = Sys.time()
  sequence_ratio = seq (min(ratio_neur_range), max(ratio_neur_range), by = step_neur_range)
  for (ratio_layers_i in sequence_ratio){
    if(!quiet) {
      cat("\ntesting with a ratio of",ratio_layers_i,"\n")
    }
    
    res_optimization_i <- find_optimal_n_neurons(spect_matrix = spect_matrix, 
                                               ref_table = ref_table,
                                               ratio_neuron_layers = ratio_layers_i,
                                               layers_range = layers_range,
                                               neurons_range = neurons_range, 
                                               step_neurons = step_neurons,
                                               quiet = T,
                                               method = method,
                                               k_fold = k_fold)
    if(res_optimization_i[1] > best_accuracy){
      best_accuracy = res_optimization_i[[1]]
      best_results = res_optimization_i
      best_ratio_layers = ratio_layers_i
      if (!quiet){
        cat("\nBest ratio for the moment:", best_ratio_layers, ", (",best_results[[3]],"neurons), with", best_accuracy*100,"% accuracy \n")
      }
    }
    
    if(!quiet){
      #time left estimation
      time_i = round(as.numeric(as.period(Sys.time()-time0, unit="sec")),0)
      perc_compl = which(ratio_layers_i == sequence_ratio)/length(sequence_ratio)
      t_tot = time_i/perc_compl
      s_left = t_tot - time_i
      m_left = as.integer(s_left / 60)
      h_left = as.integer(m_left / 60)
      s_left = as.integer(s_left %% 60)
      m_left = as.integer(m_left %% 60)
      cat("Advanced ", paste0(round(perc_compl*100,2),"%")," time left approx.:", paste0(h_left,"h:",m_left,"m:",s_left,"s"),"\n")
    }
  }
  
  best_results <- c(best_results, best_ratio_layers)
  if(!quiet){
    cat("\nDONE. BEST RESULTS:\n")
    print(best_results)
  }

  return(best_results)
}


find_optimal_nnet_dimensions_multicore <- function (cpus=3,
                                          ref_table, # reference table for the spectra
                                          spect_matrix, # matrix containing spectra
                                          seed = 222, # random seed
                                          method ="t-t",
                                          layers_range = c(1:3), # min and max number of hidden layers to test
                                          neurons_range = c(60:80), # min and max number of neurons on the first hidden layer to test
                                          step_neurons = 2, # step of increase from min number of neurons to test
                                          ratio_neur_range = c(0.1,1.0), # min and max value of ratio between hidden layers to test
                                          step_neur_range = 0.1,
                                          k_fold = 10,
                                          quiet=F) {
  if(!quiet){
    cat("looking for optiman nnet dimensions, starting: ", as.character(Sys.time()))
    cat("\nseed =",seed,", method =",method,",layers range =", layers_range,",neurons range = (", min(neurons_range),",",max(neurons_range),")\n")
    cat("step neurons =",step_neurons,",ratio neur range = (", min(ratio_neur_range),",",max(ratio_neur_range),"), step neur range =", step_neur_range)
  }
  
  
  best_accuracy = 0
  best_ratio_layers = 0
  time0 = Sys.time()
  sequence_ratio = seq (min(ratio_neur_range), max(ratio_neur_range), by = step_neur_range)
  for (ratio_layers_i in sequence_ratio){
    if(!quiet) {
      cat("\ntesting with a ratio of", ratio_layers_i, "\n")
    }
    
    res_optimization_i <- find_optimal_n_neurons_multicore(
                                                 cpus=cpus,
                                                 spect_matrix = spect_matrix, 
                                                 ref_table = ref_table,
                                                 ratio_neuron_layers = ratio_layers_i,
                                                 layers_range = layers_range,
                                                 neurons_range = neurons_range, 
                                                 step_neurons = step_neurons,
                                                 quiet = T,
                                                 method = method,
                                                 k_fold = k_fold)
    if(res_optimization_i[1] > best_accuracy){
      best_accuracy = res_optimization_i[[1]]
      best_results = res_optimization_i
      best_ratio_layers = ratio_layers_i
      if (!quiet){
        cat("\nBest ratio for the moment:", best_ratio_layers, ", (",best_results[[3]],"neurons), with", best_accuracy*100,"% accuracy \n")
      }
    }
    
    if(!quiet){
      #time left estimation
      time_i = round(as.numeric(as.period(Sys.time()-time0, unit="sec")),0)
      perc_compl = which(ratio_layers_i == sequence_ratio)/length(sequence_ratio)
      t_tot = time_i/perc_compl
      s_left = t_tot - time_i
      m_left = as.integer(s_left / 60)
      h_left = as.integer(m_left / 60)
      s_left = as.integer(s_left %% 60)
      m_left = as.integer(m_left %% 60)
      cat("Advanced ", paste0(round(perc_compl*100,2),"%")," time left approx.:", paste0(h_left,"h:",m_left,"m:",s_left,"s"),"\n")
    }
  }
  
  best_results <- c(best_results, best_ratio_layers)
  if(!quiet){
    cat("\nDONE. BEST RESULTS:\n")
    print(best_results)
  }
  
  return(best_results)
}


################################################################## END

 cat("Functions for IR spectra analysis loaded: \n",lsf.str())
