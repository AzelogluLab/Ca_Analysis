# Compare the .csv files for KO and WT in the Ca assay


WT_filenames <- list.files("Desktop/Calcium/MPodoNEBL Combined Ca Imaging", pattern = "WT.*.csv", full.names = TRUE);
KO_filenames <- list.files("Desktop/Calcium/MPodoNEBL Combined Ca Imaging", pattern = "KO.*.csv", full.names = TRUE);
WT_NumFiles <- length(WT_filenames);
KO_NumFiles <- length(KO_filenames);

WT_data <- matrix(list(), nrow = 1, ncol = WT_NumFiles);
WT_cell_count_vector <- vector(length = WT_NumFiles);
WT_size_vector <- vector(length = WT_NumFiles);
for(i in 1:WT_NumFiles){
  WT_data[[1,i]] <- read.csv(WT_filenames[i], header = FALSE, stringsAsFactors = FALSE)
  WT_cell_count_vector[i] <- (dim(WT_data[[1,i]])[2] - 1)/6
  WT_size_vector[i] <- dim(WT_data[[1,i]])[1]
}
  
KO_data <- matrix(list(), nrow = 1, ncol = KO_NumFiles);
KO_cell_count_vector <- vector(length = KO_NumFiles);
KO_size_vector <- vector(length = KO_NumFiles);
for(i1 in 1:KO_NumFiles){
  KO_data[[1,i1]] <- read.csv(KO_filenames[i1], header = FALSE, stringsAsFactors = FALSE)
  KO_cell_count_vector[i1] <- (dim(KO_data[[1,i1]])[2] - 1)/6
  KO_size_vector[i1] <- dim(KO_data[[1,i1]])[1]
}

# Storage
KO_average_nuclear_intensity <- matrix(data = NA, nrow = 600, ncol = KO_NumFiles ) ;
KO_average_cytoplasmic_intensity <- matrix(data = NA, nrow = 600, ncol = KO_NumFiles) ;
WT_average_nuclear_intensity <- matrix(data = NA, nrow = 600, ncol = WT_NumFiles) ;
WT_average_cytoplasmic_intensity <- matrix(data = NA, nrow = 600, ncol = WT_NumFiles) ;
KO_normalized_average_nuclear <- matrix(data = NA, nrow = 600, ncol = KO_NumFiles) ;
dKO_norm_avg_nuc <- vector(length = 600);
KO_normalized_average_cytoplasmic <- matrix(data = NA, nrow = 600, ncol = KO_NumFiles) ;
dKO_norm_avg_cyt <- vector(length = 600);
WT_normalized_average_nuclear <- matrix(data = NA, nrow = 600, ncol = WT_NumFiles) ;
dWT_norm_avg_nuc <- vector(length = 600);
WT_normalized_average_cytoplasmic <- matrix(data = NA, nrow = 600, ncol = WT_NumFiles) ;
dWT_norm_avg_cyt <- vector(length = 600);
KO_std_dev_max_cyt <- matrix(data = NA, nrow = 600, ncol = KO_NumFiles) ;
KO_std_dev_min_cyt <- matrix(data = NA, nrow = 600, ncol = KO_NumFiles) ;
KO_std_dev_max_nuc <- matrix(data = NA, nrow = 600, ncol = KO_NumFiles) ;
KO_std_dev_min_nuc <- matrix(data = NA, nrow = 600, ncol = KO_NumFiles) ;
WT_std_dev_max_cyt <- matrix(data = NA, nrow = 600, ncol = WT_NumFiles) ;
WT_std_dev_min_cyt <- matrix(data = NA, nrow = 600, ncol = WT_NumFiles) ;
WT_std_dev_max_nuc <- matrix(data = NA, nrow = 600, ncol = WT_NumFiles) ;
WT_std_dev_min_nuc <- matrix(data = NA, nrow = 600, ncol = WT_NumFiles) ;
KO_max_val_nuclear <- vector(length = KO_NumFiles) ;
WT_max_val_nuclear <- vector(length = WT_NumFiles) ;
KO_max_val_cytoplasmic <- vector(length = KO_NumFiles) ;
WT_max_val_cytoplasmic <- vector(length = WT_NumFiles) ;
WT_init_nuclear <- vector(length = WT_NumFiles) ;
WT_max_nuclear <- vector(length = WT_NumFiles) ;
WT_final_nuclear <- vector(length = WT_NumFiles) ;
WT_init_cytoplasmic <- vector(length = WT_NumFiles) ;
WT_max_cytoplasmic <- vector(length = WT_NumFiles) ;
WT_final_cytoplasmic <- vector(length = WT_NumFiles) ;
KO_init_nuclear <- vector(length = KO_NumFiles) ;
KO_max_nuclear <- vector(length = KO_NumFiles) ;
KO_final_nuclear <- vector(length = KO_NumFiles) ;
KO_init_cytoplasmic <- vector(length = KO_NumFiles) ;
KO_max_cytoplasmic <- vector(length = KO_NumFiles) ;
KO_final_cytoplasmic <- vector(length = KO_NumFiles) ;
time <- seq(1:600) ;

# WT data cycling and manipulation

for(y in 1:WT_NumFiles){
  sizeWT <- WT_size_vector[y]
  widthWT <- WT_cell_count_vector[y]
  WT_mean_vec_nuclear <- vector(length = widthWT) ;
  WT_mean_vec_cytoplasmic <- vector(length = widthWT) ;
  WT_max_vec_nuclear <- vector(length = widthWT) ;
  WT_min_vec_nuclear <- vector(length = widthWT) ;
  WT_max_vec_cytoplasmic <- vector(length = widthWT) ;
  WT_min_vec_cytoplasmic <- vector(length = widthWT) ;
  for(g in 1:(sizeWT)){
    for(h in 1:(widthWT)){
      WT_mean_vec_nuclear[h] <- WT_data[[1,y]][g,((6*h) - 4)]
      WT_mean_vec_cytoplasmic[h] <- WT_data[[1,y]][g,((6*h) - 1)]
      WT_max_vec_nuclear[h] <- WT_data[[1,y]][g,((6*h) - 2)]
      WT_min_vec_nuclear[h] <- WT_data[[1,y]][g,((6*h) - 3)]
      WT_max_vec_cytoplasmic[h] <- WT_data[[1,y]][g,((6*h) + 1)]
      WT_min_vec_cytoplasmic[h] <- WT_data[[1,y]][g,(6*h)]
    }
    WT_average_nuclear_intensity[g,y] <- (sum(WT_mean_vec_nuclear) / widthWT)
    WT_average_cytoplasmic_intensity[g,y] <- (sum(WT_mean_vec_cytoplasmic) / widthWT)
    WT_std_dev_max_nuc[g,y] <- (sd(WT_max_vec_nuclear)/max(WT_max_vec_nuclear))
    WT_std_dev_min_nuc[g,y] <- (sd(WT_min_vec_nuclear)/max(WT_min_vec_nuclear))
    WT_std_dev_max_cyt[g,y] <- (sd(WT_max_vec_cytoplasmic)/max(WT_max_vec_cytoplasmic))
    WT_std_dev_min_cyt[g,y] <- (sd(WT_min_vec_cytoplasmic)/max(WT_min_vec_cytoplasmic))
  }
  WT_max_val_nuclear[y] <- max(WT_average_nuclear_intensity[,y]) ;
  WT_max_val_cytoplasmic[y] <- max(WT_average_cytoplasmic_intensity[,y]) ;
  if(WT_max_val_nuclear[y] > WT_max_val_cytoplasmic[y]){
    WT_normalized_average_nuclear[,y] <- (WT_average_nuclear_intensity[,y] / WT_max_val_nuclear[y]) 
    WT_normalized_average_cytoplasmic[,y] <- (WT_average_cytoplasmic_intensity[,y] / WT_max_val_nuclear[y]) 
  } else{
    WT_normalized_average_nuclear[,y] <- (WT_average_nuclear_intensity[,y] / WT_max_val_cytoplasmic[y]) 
    WT_normalized_average_cytoplasmic[,y] <- (WT_average_cytoplasmic_intensity[,y] / WT_max_val_cytoplasmic[y]) 
  }
}

WT_avg_nuclear <- .rowMeans(WT_normalized_average_nuclear, 600, WT_NumFiles, na.rm = FALSE) ;
WT_sd_max_nuc <- .rowMeans(WT_std_dev_max_nuc, 600, WT_NumFiles, na.rm = FALSE) ;
WT_sd_min_nuc <- .rowMeans(WT_std_dev_min_nuc, 600, WT_NumFiles, na.rm = FALSE) ;
WT_avg_cytoplasmic <- .rowMeans(WT_normalized_average_cytoplasmic, 600, WT_NumFiles, na.rm = FALSE) ;
WT_sd_max_cyt <- .rowMeans(WT_std_dev_max_cyt, 600, WT_NumFiles, na.rm = FALSE) ;
WT_sd_min_cyt <- .rowMeans(WT_std_dev_min_cyt, 600, WT_NumFiles, na.rm = FALSE) ;

# KO data cycling and manipulation

for(x in 1:KO_NumFiles){
  size <- KO_size_vector[x]
  width <- KO_cell_count_vector[x]
  mean_vec_nuclear <- vector(length = width) ;
  mean_vec_cytoplasmic <- vector(length = width) ;
  max_vec_nuclear <- vector(length = width) ;
  min_vec_nuclear <- vector(length = width) ;
  max_vec_cytoplasmic <- vector(length = width) ;
  min_vec_cytoplasmic <- vector(length = width) ;
  for(i in 1:(size)){
    for(j in 1:width){
      mean_vec_nuclear[j] <- KO_data[[1,x]][i,((6*j) - 4)]
      mean_vec_cytoplasmic[j] <- KO_data[[1,x]][i,((6*j) - 1)]
      max_vec_nuclear[j] <- KO_data[[1,x]][i,((6*j) - 2)]
      min_vec_nuclear[j] <- KO_data[[1,x]][i,((6*j) - 3)]
      max_vec_cytoplasmic[j] <- KO_data[[1,x]][i,((6*j) + 1)]
      min_vec_cytoplasmic[j] <- KO_data[[1,x]][i,(6*j)]
    }
    KO_average_nuclear_intensity[i,x] <- (sum(mean_vec_nuclear) / width)
    KO_average_cytoplasmic_intensity[i,x] <- (sum(mean_vec_cytoplasmic) / width)
    KO_std_dev_max_nuc[i,x] <- (sd(max_vec_nuclear)/max(max_vec_nuclear))
    KO_std_dev_min_nuc[i,x] <- (sd(min_vec_nuclear)/max(min_vec_nuclear))
    KO_std_dev_max_cyt[i,x] <- (sd(max_vec_cytoplasmic)/max(max_vec_cytoplasmic))
    KO_std_dev_min_cyt[i,x] <- (sd(min_vec_cytoplasmic)/max(min_vec_cytoplasmic))
  }
  KO_max_val_nuclear[x] <- max(KO_average_nuclear_intensity[,x]) ;
  KO_max_val_cytoplasmic[x] <- max(KO_average_cytoplasmic_intensity[,x]) ;
  if(KO_max_val_nuclear[x] > KO_max_val_cytoplasmic[x]){
    KO_normalized_average_nuclear[,x] <- (KO_average_nuclear_intensity[,x] / KO_max_val_nuclear[x]) 
    KO_normalized_average_cytoplasmic[,x] <- (KO_average_cytoplasmic_intensity[,x] / KO_max_val_nuclear[x]) 
  } else{
    KO_normalized_average_nuclear[,x] <- (KO_average_nuclear_intensity[,x] / KO_max_val_cytoplasmic[x]) 
    KO_normalized_average_cytoplasmic[,x] <- (KO_average_cytoplasmic_intensity[,x] / KO_max_val_cytoplasmic[x]) 
  }
}

KO_avg_nuclear <- .rowMeans(KO_normalized_average_nuclear, 600, KO_NumFiles, na.rm = FALSE) ;
KO_sd_max_nuc <- .rowMeans(KO_std_dev_max_nuc, 600, KO_NumFiles, na.rm = FALSE) ;
KO_sd_min_nuc <- .rowMeans(KO_std_dev_min_nuc, 600, KO_NumFiles, na.rm = FALSE) ;
KO_avg_cytoplasmic <- .rowMeans(KO_normalized_average_cytoplasmic, 600, KO_NumFiles, na.rm = FALSE) ;
KO_sd_max_cyt <- .rowMeans(KO_std_dev_max_cyt, 600, KO_NumFiles, na.rm = FALSE) ;
KO_sd_min_cyt <- .rowMeans(KO_std_dev_min_cyt, 600, KO_NumFiles, na.rm = FALSE) ;

# Calculation of derivatives (aggregate and individual)

for(s in 2:600){
  dWT_norm_avg_nuc[s] <- (WT_avg_nuclear[s] - WT_avg_nuclear[s-1])
  dWT_norm_avg_cyt[s] <- (WT_avg_cytoplasmic[s] - WT_avg_cytoplasmic[s-1])
  dKO_norm_avg_nuc[s] <- (KO_avg_nuclear[s] - KO_avg_nuclear[s-1])
  dKO_norm_avg_cyt[s] <- (KO_avg_cytoplasmic[s] - KO_avg_cytoplasmic[s-1])
}

max_aggregate_deriv <- vector(length = 4)
max_aggregate_deriv[1] <- max(dWT_norm_avg_cyt)
max_aggregate_deriv[2] <- max(dWT_norm_avg_nuc)
max_aggregate_deriv[3] <- max(dKO_norm_avg_cyt)
max_aggregate_deriv[4] <- max(dKO_norm_avg_nuc)

dWT_cyt_trials <- matrix(data = NA, nrow = 599, ncol = WT_NumFiles)
dWT_nuc_trials <- matrix(data = NA, nrow = 599, ncol = WT_NumFiles)
dKO_cyt_trials <- matrix(data = NA, nrow = 599, ncol = KO_NumFiles)
dKO_nuc_trials <- matrix(data = NA, nrow = 599, ncol = KO_NumFiles)
dWT_cyt_max <- vector(length = WT_NumFiles)
dWT_nuc_max <- vector(length = WT_NumFiles)
dKO_cyt_max <- vector(length = WT_NumFiles)
dKO_nuc_max <- vector(length = WT_NumFiles)

for(i6 in 1:WT_NumFiles){
  for(i7 in 2:600){
    dWT_cyt_trials[[i7-1,i6]] <- (WT_normalized_average_cytoplasmic[[i7,i6]] - WT_normalized_average_cytoplasmic[[i7-1,i6]])
    dWT_nuc_trials[[i7-1,i6]] <- (WT_normalized_average_nuclear[[i7,i6]] - WT_normalized_average_nuclear[[i7-1,i6]])
  }
  dWT_cyt_max[i6] <- max(dWT_cyt_trials[,i6])
  dWT_nuc_max[i6] <- max(dWT_nuc_trials[,i6])
}

for(i8 in 1:KO_NumFiles){
  for(i9 in 2:600){
    dKO_cyt_trials[[i9-1,i8]] <- (KO_normalized_average_cytoplasmic[[i9,i8]] - KO_normalized_average_cytoplasmic[[i9-1,i8]])
    dKO_nuc_trials[[i9-1,i8]] <- (KO_normalized_average_nuclear[[i9,i8]] - KO_normalized_average_nuclear[[i9-1,i8]])
  }
  dKO_cyt_max[i8] <- max(dKO_cyt_trials[,i8])
  dKO_nuc_max[i8] <- max(dKO_nuc_trials[,i8])
}
dKO_cyt_max[14:16] <- "NA"
dKO_nuc_max[14:16] <- "NA"
# Visual Representation of Time-Series Data + Error Bars

par(mfrow = c(2,2))
transparent_khaki <- adjustcolor("khaki", alpha.f = 0.5) ;
transparent_blue <- adjustcolor("skyblue", alpha.f = 0.3) ;
transparent_green <- adjustcolor("light green", alpha.f = 0.5) ;
transparent_cyan <- adjustcolor("cyan", alpha.f = 0.3) ;

plot(time, KO_avg_cytoplasmic, col = "green", xlim = range(time), ylim = range(c(-0.3,1.3)), ylab = "Mean Normalized Intensity", main = "Response to
     Calcium Stimulus (KO)", xlab = "Time", type = "l")
arrows(time, KO_avg_cytoplasmic + KO_sd_max_cyt, time, KO_avg_cytoplasmic - KO_sd_min_cyt, length = 0.001, angle = 90, code = 3, col = transparent_green)
lines(time[order(time)], KO_avg_nuclear[order(time)], col = "blue", xlim = range(time), ylim = range(c(-0.3,1.3)))
arrows(time, KO_avg_nuclear + KO_sd_max_nuc, time, KO_avg_nuclear - KO_sd_min_nuc, length = 0.001, angle = 90, code = 3, col = transparent_cyan)
legend(300, 1.3, legend = c("Nucleus", "Cytoplasm"), col = c("blue", "green"), lty = c(1,1), cex = 0.6)

plot(time, WT_avg_cytoplasmic, col = "green", xlim = range(time), ylim = range(c(-0.3,1.3)), ylab = "Mean Normalized Intensity", main = "Response to
     Calcium Stimulus (WT)", xlab = "Time", type = "l")
arrows(time, WT_avg_cytoplasmic + WT_sd_max_cyt, time, WT_avg_cytoplasmic - WT_sd_min_cyt, length = 0.001, angle = 90, code = 3, col = transparent_green)
lines(time[order(time)], WT_avg_nuclear[order(time)], col = "blue", xlim = range(time), ylim = range(c(-0.3,1.3)))
arrows(time, WT_avg_nuclear + WT_sd_max_nuc, time, WT_avg_nuclear - WT_sd_min_nuc, length = 0.001, angle = 90, code = 3, col = transparent_cyan)
legend(300, 1.3, legend = c("Nucleus", "Cytoplasm"), col = c("blue", "green"), lty = c(1,1), cex = 0.6)

plot(time, KO_avg_cytoplasmic, col = "orange", xlim = range(time), ylim = range(c(-0.3,1.3)), ylab = "Mean Normalized Intensity", main = "Cytoplasmic Calcium 
     Response", xlab = "Time", type = "l")
arrows(time, KO_avg_cytoplasmic + KO_sd_max_cyt, time, KO_avg_cytoplasmic - KO_sd_min_cyt, length = 0.001, angle = 90, code = 3, col = transparent_khaki)
lines(time[order(time)], WT_avg_cytoplasmic[order(time)], col = "blue", xlim = range(time), ylim = range(c(-0.3,1.3)))
arrows(time, WT_avg_cytoplasmic + WT_sd_max_cyt, time, WT_avg_cytoplasmic - WT_sd_min_cyt, length = 0.001, angle = 90, code = 3, col = transparent_blue)
legend(300, 1.3, legend = c("WT", "KO"), col = c("blue", "khaki"), lty = c(1,1), cex = 0.6)

plot(time, KO_avg_nuclear, col = "orange", xlim = range(time), ylim = range(c(-0.3,1.3)), ylab = "Mean Normalized Intensity", main = "Nuclear Calcium 
     Response", xlab = "Time", type = "l")
arrows(time, KO_avg_nuclear + KO_sd_max_nuc, time, KO_avg_nuclear - KO_sd_min_nuc, length = 0.001, angle = 90, code = 3, col = transparent_khaki)
lines(time[order(time)], WT_avg_nuclear[order(time)], col = "blue", xlim = range(time), ylim = range(c(-0.3,1.3)))
arrows(time, WT_avg_nuclear + WT_sd_max_nuc, time, WT_avg_nuclear - WT_sd_min_nuc, length = 0.001, angle = 90, code = 3, col = transparent_blue)
legend(300, 1.3, legend = c("WT", "KO"), col = c("blue", "khaki"), lty = c(1,1), cex = 0.6)

par(mfrow = c(2,1))
plot(time, dKO_norm_avg_nuc, ylab = "d[Intensity]/dt", xlab = "Time", main = "Nuclear Derivative", type = "l", col = "khaki")
lines(time[order(time)], dWT_norm_avg_nuc[order(time)], col = "blue")
legend(400, max(dKO_norm_avg_nuc), legend = c("WT", "KO"), col = c("blue", "khaki"), lty = c(1,1), cex = 0.6)

plot(time, dKO_norm_avg_cyt, ylab = "d[Intensity]/dt", xlab = "Time", main = "Cytoplasmic Derivative", type = "l", col = "khaki")
lines(time[order(time)], dWT_norm_avg_cyt[order(time)], col = "blue")
legend(400, max(dKO_norm_avg_cyt), legend = c("WT", "KO"), col = c("blue", "khaki"), lty = c(1,1), cex = 0.6)

# Normalized values adjusted to reduce difference in initial values and visual presentation of adjusted time-series data

par(mfrow = c(1,1))

Adj_KO_avg_cytoplasmic <- vector(length = 600)
Adj_KO_avg_nuclear <- vector(length = 600)
Adj_WT_avg_nuclear <- vector(length = 600)
Adj_WT_avg_cytoplasmic <- vector(length = 600)

if(KO_avg_nuclear[1] > WT_avg_nuclear[1]){
  for(i2 in 1:600){
    Diff_nuc <- KO_avg_nuclear[1] - WT_avg_nuclear[1]
    Adj_KO_avg_nuclear[i2] <- KO_avg_nuclear[i2] - Diff_nuc
  }
  plot(time, Adj_KO_avg_nuclear, col = rgb(251,209,120,255, names = NULL, maxColorValue = 255), xlim = range(time), ylim = range(c(-0.3,1.3)), ylab = "Mean Normalized Intensity", main = "Nuclear Calcium 
       Response", xlab = "Time", type = "l")
  arrows(time, Adj_KO_avg_nuclear + KO_sd_max_nuc, time, Adj_KO_avg_nuclear - KO_sd_min_nuc, length = 0.001, angle = 90, code = 3, col = transparent_khaki)
  lines(time[order(time)], WT_avg_nuclear[order(time)], col = rgb(0,83,144,255, names = NULL, maxColorValue = 255), xlim = range(time), ylim = range(c(-0.3,1.3)))
  arrows(time, WT_avg_nuclear + WT_sd_max_nuc, time, WT_avg_nuclear - WT_sd_min_nuc, length = 0.001, angle = 90, code = 3, col = transparent_blue)
  legend(300, 1.3, legend = c("WT", "KO"), col = c(rgb(0,83,144,255, names = NULL, maxColorValue = 255), rgb(251,209,120,255, names = NULL, maxColorValue = 255)), lty = c(1,1), cex = 0.6)
} else{
  for(i2 in 1:600){
    Diff_nuc <- WT_avg_nuclear[1] - KO_avg_nuclear[1]
    Adj_WT_avg_nuclear[i2] <- WT_avg_nuclear[i2] - Diff_nuc
  }
  plot(time, KO_avg_nuclear, col = rgb(251,209,120,255, names = NULL, maxColorValue = 255), xlim = range(time), ylim = range(c(-0.3,1.3)), ylab = "Mean Normalized Intensity", main = "Nuclear Calcium 
       Response", xlab = "Time", type = "l")
  arrows(time, KO_avg_nuclear + KO_sd_max_nuc, time, KO_avg_nuclear - KO_sd_min_nuc, length = 0.001, angle = 90, code = 3, col = transparent_khaki)
  lines(time[order(time)], Adj_WT_avg_nuclear[order(time)], col = rgb(0,83,144,255, names = NULL, maxColorValue = 255), xlim = range(time), ylim = range(c(-0.3,1.3)))
  arrows(time, Adj_WT_avg_nuclear + WT_sd_max_nuc, time, Adj_WT_avg_nuclear - WT_sd_min_nuc, length = 0.001, angle = 90, code = 3, col = transparent_blue)
  legend(300, 1.3, legend = c("WT", "KO"), col = c(rgb(0,83,144,255, names = NULL, maxColorValue = 255), rgb(251,209,120,255, names = NULL, maxColorValue = 255)), lty = c(1,1), cex = 0.6)
}

if(KO_avg_cytoplasmic[1] > WT_avg_cytoplasmic[1]){
  for(i3 in 1:600){
    Diff_cyt <- KO_avg_cytoplasmic[1] - WT_avg_cytoplasmic[1]
    Adj_KO_avg_cytoplasmic[i3] <- KO_avg_cytoplasmic[i3] - Diff_cyt
  }
  plot(time, Adj_KO_avg_cytoplasmic, col = rgb(251,209,120,255, names = NULL, maxColorValue = 255), xlim = range(time), ylim = range(c(-0.3,1.3)), ylab = "Mean Normalized Intensity", main = "Cytoplasmic Calcium 
       Response", xlab = "Time", type = "l")
  arrows(time, Adj_KO_avg_cytoplasmic + KO_sd_max_cyt, time, Adj_KO_avg_cytoplasmic - KO_sd_min_cyt, length = 0.001, angle = 90, code = 3, col = transparent_khaki)
  lines(time[order(time)], WT_avg_cytoplasmic[order(time)], col = rgb(0,83,144,255, names = NULL, maxColorValue = 255), xlim = range(time), ylim = range(c(-0.3,1.3)))
  arrows(time, WT_avg_cytoplasmic + WT_sd_max_cyt, time, WT_avg_cytoplasmic - WT_sd_min_cyt, length = 0.001, angle = 90, code = 3, col = transparent_blue)
  legend(300, 1.3, legend = c("WT", "KO"), col = c(rgb(0,83,144,255, names = NULL, maxColorValue = 255), rgb(251,209,120,255, names = NULL, maxColorValue = 255)), lty = c(1,1), cex = 0.6)
} else{
  for(i3 in 1:600){
    Diff_cyt <- WT_avg_cytoplasmic[1] - KO_avg_cytoplasmic[1]
    Adj_WT_avg_cytoplasmic[i3] <- WT_avg_cytoplasmic[i3] - Diff_cyt
  }
  plot(time, KO_avg_cytoplasmic, col = rgb(251,209,120,255, names = NULL, maxColorValue = 255), xlim = range(time), ylim = range(c(-0.3,1.3)), ylab = "Mean Normalized Intensity", main = "Cytoplasmic Calcium 
       Response", xlab = "Time", type = "l")
  arrows(time, KO_avg_cytoplasmic + KO_sd_max_cyt, time, KO_avg_cytoplasmic - KO_sd_min_cyt, length = 0.001, angle = 90, code = 3, col = transparent_khaki)
  lines(time[order(time)], Adj_WT_avg_cytoplasmic[order(time)], col = rgb(0,83,144,255, names = NULL, maxColorValue = 255), xlim = range(time), ylim = range(c(-0.3,1.3)))
  arrows(time, Adj_WT_avg_cytoplasmic + WT_sd_max_cyt, time, Adj_WT_avg_cytoplasmic - WT_sd_min_cyt, length = 0.001, angle = 90, code = 3, col = transparent_blue)
  legend(300, 1.3, legend = c("WT", "KO"), col = c(rgb(0,83,144,255, names = NULL, maxColorValue = 255), rgb(251,209,120,255, names = NULL, maxColorValue = 255)), lty = c(1,1), cex = 0.6)
}


# Curve Fitting
library(nlme)
# KO Data

max_KO_cyt <- max(KO_avg_cytoplasmic);
max_KO_nuc <- max(KO_avg_nuclear);
KO_cyt_ind <- max_KO_cyt == KO_avg_cytoplasmic;
KO_nuc_ind <- max_KO_nuc == KO_avg_nuclear;
KO_cyt_index <- match(KO_avg_cytoplasmic[KO_cyt_ind],KO_avg_cytoplasmic)
KO_nuc_index <- match(KO_avg_nuclear[KO_nuc_ind],KO_avg_nuclear)


KO_fitting_data <- matrix(list(), nrow = 2, ncol = KO_NumFiles)

for(z in 1:KO_NumFiles){
  sizeKO <- KO_size_vector[z]
  widthKO <- KO_cell_count_vector[z]
  KO_nuc_mat <- matrix(data = NA, nrow = sizeKO, ncol = widthKO)
  KO_cyt_mat <- matrix(data = NA, nrow = sizeKO, ncol = widthKO)
  for(f in 1:widthKO){
    KO_nuc_mat[,f] <- KO_data[[1,z]][,((6*f) - 4)]
    KO_cyt_mat[,f] <- KO_data[[1,z]][,((6*f) - 1)]
    KO_nuc_max <- max(KO_nuc_mat[,f])
    KO_cyt_max <- max(KO_cyt_mat[,f])
    KO_nuc_mat[,f] <- KO_nuc_mat[,f]/KO_nuc_max
    KO_cyt_mat[,f] <- KO_cyt_mat[,f]/KO_cyt_max
  }
  KO_fitting_data[[1,z]] <- KO_nuc_mat
  KO_fitting_data[[2,z]] <- KO_cyt_mat
}
KO_fit_df <- as.data.frame(KO_fitting_data) ;
KO_fit_coef <- matrix(list(), nrow = 2, ncol = KO_NumFiles)
KO_half_life_data <- matrix(list(), nrow = 2, ncol = KO_NumFiles)

for(z1 in 1:KO_NumFiles){
  sizeKO <- KO_size_vector[z1]
  widthKO <- KO_cell_count_vector[z1]
  KO_fit_vector_nuc <- vector(length = widthKO)
  KO_fit_vector_cyt <- vector(length = widthKO)
  KO_half_life_nuc <- vector(length = widthKO)
  KO_half_life_cyt <- vector(length = widthKO)
  for(f1 in 1:widthKO){
    ds_nuc <- data.frame(x = time[1:(250 - KO_nuc_index + 1)], y = KO_fit_df[[1,z1]][KO_nuc_index:250,f1])
    ds_cyt <- data.frame(x = time[1:(250 - KO_cyt_index + 1)], y = KO_fit_df[[2,z1]][KO_cyt_index:250,f1])
    nuc_fit <- gnls(y ~ I(a*exp(b*x)), data = ds_nuc, start = list(a=1, b=0))
    cyt_fit <- gnls(y ~ I(a*exp(b*x)), data = ds_cyt, start = list(a=1, b=0))
    KO_fit_vector_nuc[f1] <- coef(nuc_fit)[2]
    KO_fit_vector_cyt[f1] <- coef(cyt_fit)[2]
    KO_half_life_nuc[f1] <- log(0.5) / coef(nuc_fit)[2]
    KO_half_life_cyt[f1] <- log(0.5) / coef(cyt_fit)[2]
  }
  KO_fit_coef[[1,z1]] <- KO_fit_vector_nuc
  KO_fit_coef[[2,z1]] <- KO_fit_vector_cyt
  KO_half_life_data[[1,z1]] <- KO_half_life_nuc
  KO_half_life_data[[2,z1]] <- KO_half_life_cyt
}

# WT Data

max_WT_cyt <- max(WT_avg_cytoplasmic);
max_WT_nuc <- max(WT_avg_nuclear);
WT_cyt_ind <- max_WT_cyt == WT_avg_cytoplasmic;
WT_nuc_ind <- max_WT_nuc == WT_avg_nuclear;
WT_cyt_index <- match(WT_avg_cytoplasmic[WT_cyt_ind],WT_avg_cytoplasmic)
WT_nuc_index <- match(WT_avg_nuclear[WT_nuc_ind],WT_avg_nuclear)

WT_fitting_data <- matrix(list(), nrow = 2, ncol = WT_NumFiles)

for(z in 1:WT_NumFiles){
  sizeWT <- WT_size_vector[z]
  widthWT <- WT_cell_count_vector[z]
  WT_nuc_mat <- matrix(data = NA, nrow = sizeWT, ncol = widthWT)
  WT_cyt_mat <- matrix(data = NA, nrow = sizeWT, ncol = widthWT)
  for(f in 1:widthWT){
    WT_nuc_mat[,f] <- WT_data[[1,z]][,((6*f) - 4)]
    WT_cyt_mat[,f] <- WT_data[[1,z]][,((6*f) - 1)]
    WT_nuc_max <- max(WT_nuc_mat[,f])
    WT_cyt_max <- max(WT_cyt_mat[,f])
    WT_nuc_mat[,f] <- WT_nuc_mat[,f]/WT_nuc_max
    WT_cyt_mat[,f] <- WT_cyt_mat[,f]/WT_cyt_max
  }
  WT_fitting_data[[1,z]] <- WT_nuc_mat
  WT_fitting_data[[2,z]] <- WT_cyt_mat
}
WT_fit_df <- as.data.frame(WT_fitting_data) ;
WT_fit_coef <- matrix(list(), nrow = 2, ncol = WT_NumFiles)
WT_half_life_data <- matrix(list(), nrow = 2, ncol = WT_NumFiles)

for(z1 in 1:WT_NumFiles){
  sizeWT <- WT_size_vector[z1]
  widthWT <- WT_cell_count_vector[z1]
  WT_fit_vector_nuc <- vector(length = widthWT)
  WT_fit_vector_cyt <- vector(length = widthWT)
  WT_half_life_nuc <- vector(length = widthWT)
  WT_half_life_cyt <- vector(length = widthWT)
  for(f1 in 1:widthWT){
    ds_nuc <- data.frame(x = time[1:(250 - WT_nuc_index + 1)], y = WT_fit_df[[1,z1]][WT_nuc_index:250,f1])
    ds_cyt <- data.frame(x = time[1:(250 - WT_cyt_index + 1)], y = WT_fit_df[[2,z1]][WT_cyt_index:250,f1])
    nuc_fit <- gnls(y ~ I(a*exp(b*x)), data = ds_nuc, start = list(a=1, b=0))
    cyt_fit <- gnls(y ~ I(a*exp(b*x)), data = ds_cyt, start = list(a=1, b=0))
    WT_fit_vector_nuc[f1] <- coef(nuc_fit)[2]
    WT_fit_vector_cyt[f1] <- coef(cyt_fit)[2]
    WT_half_life_nuc[f1] <- log(0.5) / coef(nuc_fit)[2]
    WT_half_life_cyt[f1] <- log(0.5) / coef(cyt_fit)[2]
  }
  WT_fit_coef[[1,z1]] <- WT_fit_vector_nuc
  WT_fit_coef[[2,z1]] <- WT_fit_vector_cyt
  WT_half_life_data[[1,z1]] <- WT_half_life_nuc
  WT_half_life_data[[2,z1]] <- WT_half_life_cyt
}

ds_KO_nuc <- data.frame(x1 = time[KO_nuc_index:250], y1 = KO_avg_nuclear[KO_nuc_index:250]) ;
ds_KO_cyt <- data.frame(x3 = time[KO_cyt_index:250], y3 = KO_avg_cytoplasmic[KO_cyt_index:250])
ds_WT_nuc <- data.frame(x2 = time[WT_nuc_index:250], y2 = WT_avg_nuclear[WT_nuc_index:250]) ;
ds_WT_cyt <- data.frame(x4 = time[WT_cyt_index:250], y4 = WT_avg_cytoplasmic[WT_cyt_index:250])
fit1 <- gnls(y1 ~ I(a*exp(b*x1)), data = ds_KO_nuc, start = list(a=1,b=0)) ;
fit3 <- gnls(y3 ~ I(a*exp(b*x3)), data = ds_KO_cyt, start = list(a=1,b=0)) ;
fit2 <- gnls(y2 ~ I(a*exp(b*x2)), data = ds_WT_nuc, start = list(a=1,b=0)) ;
fit4 <- gnls(y4 ~ I(a*exp(b*x4)), data = ds_WT_cyt, start = list(a=1,b=0)) ;

half_lives <- data.frame(KO_nuclear = (log(0.5)/coef(fit1)[2]), KO_cytoplasmic = (log(0.5)/coef(fit3)[2]), WT_nuclear = (log(0.5)/coef(fit2)[2]), 
                         WT_cytoplasmic = (log(0.5)/coef(fit4)[2]), row.names = "Half-Life");


par(mfrow=c(2,2))

plot(ds_WT_nuc$x2, ds_WT_nuc$y2, main = "WT Nuclear Fit", ylim = range(c(0,1)), xlab = "Time", ylab = "Mean Normalized Intensity")
points(ds_WT_nuc$x2, (coef(fit2)[1] * exp(coef(fit2)[2] * ds_WT_nuc$x2)), col = "red")
legend(200, 1, legend = c("Data", "Fit"), col = c("black", "red"), lty = c(1,1), cex = 0.6)

plot(ds_WT_cyt$x4, ds_WT_cyt$y4, main = "WT Cytoplasmic Fit", ylim = range(c(0,1)), xlab = "Time", ylab = "Mean Normalized Intensity")
points(ds_WT_cyt$x4, (coef(fit4)[1] * exp(coef(fit4)[2] * ds_WT_cyt$x4)), col = "red")
legend(200, 1, legend = c("Data", "Fit"), col = c("black", "red"), lty = c(1,1), cex = 0.6)

plot(ds_KO_nuc$x1, ds_KO_nuc$y1, main = "KO Nuclear Fit", ylim = range(c(0,1)), xlab = "Time", ylab = "Mean Normalized Intensity")
points(ds_KO_nuc$x1, (coef(fit1)[1] * exp(coef(fit1)[2] * ds_KO_nuc$x1)), col = "red")
legend(200, 1, legend = c("Data", "Fit"), col = c("black", "red"), lty = c(1,1), cex = 0.6)

plot(ds_KO_cyt$x3, ds_KO_cyt$y3, main = "KO Cytoplasmic Fit", ylim = range(c(0,1)), xlab = "Time", ylab = "Mean Normalized Intensity")
points(ds_KO_cyt$x3, (coef(fit3)[1] * exp(coef(fit3)[2] * ds_KO_cyt$x3)), col = "red")
legend(200, 1, legend = c("Data", "Fit"), col = c("black", "red"), lty = c(1,1), cex = 0.6)

# Statistics
# ANOVA

library(nlme)

avg_data <- data.frame(time, WT_avg_cytoplasmic, KO_avg_cytoplasmic, WT_avg_nuclear, KO_avg_nuclear) ;
stacked <- data.frame(avg_data[1], stack(avg_data[2:5])) ;

library(ggplot2)

qplot(ind, values, data = stacked, geom = "boxplot", xlab = "Group", ylab = "Mean Normalized Intensity", main = "Comparison of Mean Normalized Intensity")
lme_intensity <- lme(values ~ ind, random = ~1 | time/ind, data = stacked, method = "ML") ;
summary(lme_intensity)
ANOVA <- anova(lme_intensity)
print(ANOVA)

# Post-Hoc Testing: Tukey

library(multcomp)

posthoc <- glht(lme_intensity, linfct = mcp(ind = "Tukey"))
summary(posthoc)

print(half_lives)

# Generate Data Frame and Write to CSV


for(i4 in 1:WT_NumFiles){
  WT_init_cytoplasmic[i4] <- WT_normalized_average_cytoplasmic[1,i4] - Diff_cyt
  WT_max_cytoplasmic[i4] <- max(WT_normalized_average_cytoplasmic[,i4]) - Diff_cyt
  WT_final_cytoplasmic[i4] <- WT_normalized_average_cytoplasmic[600,i4] - Diff_cyt
  WT_init_nuclear[i4] <- WT_normalized_average_nuclear[1,i4] - Diff_nuc
  WT_max_nuclear[i4] <- max(WT_normalized_average_nuclear[,i4]) - Diff_nuc
  WT_final_nuclear[i4] <- WT_normalized_average_nuclear[600,i4] - Diff_nuc
}


for(i5 in 1:KO_NumFiles){
  KO_init_cytoplasmic[i5] <- KO_normalized_average_cytoplasmic[1,i5]
  KO_max_cytoplasmic[i5] <- max(KO_normalized_average_cytoplasmic[,i5])
  KO_final_cytoplasmic[i5] <- KO_normalized_average_cytoplasmic[600,i5]
  KO_init_nuclear[i5] <- KO_normalized_average_nuclear[1,i5]
  KO_max_nuclear[i5] <- max(KO_normalized_average_nuclear[,i5])
  KO_final_nuclear[i5] <- KO_normalized_average_nuclear[600,i5]
}

adjustment_cyt <- mean(WT_init_cytoplasmic) - mean(KO_init_cytoplasmic)
adjustment_nuc <- mean(WT_init_nuclear) - mean(KO_init_nuclear)

WT_init_max_final_data <- data.frame(WT_init_cytoplasmic, WT_max_cytoplasmic, WT_final_cytoplasmic, WT_init_nuclear,
                                     WT_max_nuclear, WT_final_nuclear)

KO_init_max_final_data <- data.frame(KO_init_cytoplasmic, KO_max_cytoplasmic, KO_final_cytoplasmic, KO_init_nuclear,
                                     KO_max_nuclear, KO_final_nuclear)

Relevant_Vals <- data.frame(KO_avg_cytoplasmic, KO_avg_nuclear, Adj_WT_avg_cytoplasmic, Adj_WT_avg_nuclear, 
                            KO_sd_max_cyt, KO_sd_min_cyt, KO_sd_max_nuc, KO_sd_min_nuc, WT_sd_max_cyt,
                            WT_sd_min_cyt, WT_sd_max_nuc, WT_sd_min_nuc, dKO_norm_avg_cyt, dKO_norm_avg_nuc, 
                            dWT_norm_avg_cyt, dWT_norm_avg_nuc)

r_names <- c("WT_Cyt Avg. Deriv","WT_Nuc Avg. Deriv","KO_Cyt Avg. Deriv","KO_Nuc Avg. Deriv")
Agg_derivs <- data.frame(max_aggregate_deriv,row.names = r_names)

WT_Trial_derivs <- data.frame(dWT_cyt_max,dWT_nuc_max)

KO_Trial_derivs <- data.frame(dKO_cyt_max,dKO_nuc_max)

Overall_Trial_derivs <- data.frame(dWT_cyt_max,dWT_nuc_max,dKO_cyt_max,dKO_nuc_max,check.rows = FALSE)

#write.csv(Relevant_Vals, file = "Relevant_Values.csv")

write.csv(WT_init_max_final_data, file = "WT_Initial_Max_Final.csv")
write.csv(KO_init_max_final_data, file = "KO_Initial_Max_Final.csv")
write.csv(Agg_derivs, file = "Max_Average_Derivative_Data.csv")
write.csv(Overall_Trial_derivs, file = "Max_Per_Trial_Derivative_Data.csv")
