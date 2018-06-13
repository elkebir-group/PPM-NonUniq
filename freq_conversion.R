# read in frequency output and convert to .data for spruce to enumerate

# add f average column
# add mu: (1,0,0) (1,1,1)

library(readr)

rewrite_tree2freq_for_downsample <- function(filename, coverage = 10000, input_folder = './tree2freqs_output/', output_folder = './downsample_input/'){
  output_file <- paste(output_folder, filename, sep = '')
  file <- read_lines(paste(input_folder, filename, sep = ''))
  if (length(file) < 2){return()}
  # file <- file[-c(1)]
  
  table <- file[c(5:length(file))]
  write_lines(table, 'temp.txt')
  mtx <- read_delim("temp.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  colnames(mtx) <- c('sample_index','sample_label','anatomical_site_index','anatomical_site_label','character_index','character_label','ref','var')
  
  mtx[,7] <- floor(coverage*(1-mtx[,8]))
  mtx[,8] <- coverage-mtx[,7]
  
  write_tsv(mtx, 'temp.txt')
  new_file <- read_lines('temp.txt')
  new_file[1] <- paste('#', new_file[1], sep = '')
  new_file <- c(file[1:3], new_file)
  # new_file <- new_file[-c(3)] #paste('#', new_file[3], sep = ' ')
  write_lines(new_file, output_file)
  return(mtx)
}

rewrite_tree2freq_for_downsample_wrapper <- function(M = c(0.1,1,10), S = c(0:30)){
  for (m in M){
    for (s in S){
      file <- paste('freqs_M', m, '_S', s, '.txt', sep = '')
      # print(file)
      rewrite_tree2freq_for_downsample(filename = file)
    }
  }
}

rewrite_downsample_for_spruce <- function(freqs_filename, downsample_filename, coverage = 10000, downsample_input_folder = './downsample_output/', freqs_input_folder = './reads2freqs_output/', output_folder = './spruce_input/'){
  output_file <- paste(output_folder, freqs_filename, sep = '')
  
  freqs_file <- read_lines(paste(freqs_input_folder, freqs_filename, sep = ''))
  if (length(freqs_file) < 2){return()}
  # freqs_file <- freqs_file[-c(1)]
  
  downsample_file <- read_lines(paste(downsample_input_folder, downsample_filename, sep = ''))
  if (length(downsample_file) < 2){return()}
  # downsample_file <- downsample_file[-c(1)]
  
  table <- downsample_file[c(5:length(downsample_file))]
  write_lines(table, 'temp.txt')
  downsample_mtx <- read_delim("temp.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  colnames(downsample_mtx) <- c('sample_index','sample_label','anatomical_site_index','anatomical_site_label','character_index','character_label','ref','var')
  f_avg <- downsample_mtx[,8]/(downsample_mtx[,8]+downsample_mtx[,7])
  
  table <- freqs_file[c(5:length(freqs_file))]
  write_lines(table, 'temp.txt')
  freqs_mtx <- read_delim("temp.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  colnames(freqs_mtx) <- c('sample_index','sample_label','anatomical_site_index','anatomical_site_label','character_index','character_label','f-','f+')
  sample_idx  <- freqs_mtx[,1]
  sample_label  <- freqs_mtx[,2]
  snv_idx   <- freqs_mtx[,5]
  snv_label <- freqs_mtx[,6]
  f_minus <- freqs_mtx[,7]
  f_plus <- freqs_mtx[,8]
  
  output_mtx <- freqs_mtx
  
  output_mtx[,1] <- sample_idx
  output_mtx[,2] <- sample_label
  output_mtx[,3] <- snv_idx
  output_mtx[,4] <- snv_idx
  output_mtx[,5] <- f_minus
  output_mtx[,6] <- f_avg
  output_mtx[,7] <- f_plus
  
  # set x,y,mu
  output_mtx[,8] <- as.integer(1)
  output_mtx[,9] <- as.integer(1)
  output_mtx[,10] <- as.integer(1)
  
  colnames(output_mtx) <- c('sample_index','sample_label','character_index','character_label','f-','f','f+','x','y','mu')
  
  write_delim(output_mtx, 'temp.txt', delim = " ")
  new_file <- read_lines('temp.txt')
  new_file[1] <- paste('#', new_file[1], sep = '')
  # print(file[1:2])
  new_file <- c(freqs_file[1:2], new_file)
  m <- strsplit(new_file[1], split = '#')[[1]][1]
  n <- strsplit(new_file[2], split = '#')[[1]][1]
  
  new_file[1] <- paste(m, '#m', sep = '')
  new_file[2] <- paste(n, '#n', sep = '')
  new_file <- new_file[-c(3)] #paste('#', new_file[3], sep = ' ')
  write_lines(new_file, output_file)
  # return(mtx)
}

rewrite_downsample_for_spruce_wrapper <- function(C = c(-1,50,100,1000), P = c(1), K = c(2,5,10), M = c(0.1,1,10), 
                                                  S = c(0:30)){
  for (c in C){
    for (p in P){
      for (k in K){
        for (m in M){
          for (s in S){
            downsample_file <- paste('downsampled_C', c, '_P', p, '_k', k, '_M', m, '_S', s, '.txt', sep = '')
            freqs_file <- paste('freq_C', c, '_P', p, '_k', k, '_M', m, '_S', s, '.txt', sep = '')
            rewrite_downsample_for_spruce(downsample_filename = downsample_file, freqs_filename = freqs_file)
          }
        }
      }
    }
  }
  
}

num_snv_analysis <- function(M = c(0.1,1,10), S = c(0:30)){
  n_snv <- matrix(data = NA, ncol = length(M), nrow = length(S))
  colnames(n_snv) <- M
  rownames(n_snv) <- S
  
  # plot.new()
  plot(0,0, xlim = c(0, 12), ylim = c(0, 1000))
  # plot(1,lwd=0,axes=c(0, 10, 0, 1000),xlab="",ylab="",...)
  
  for (m in 1:length(M)){
    for (s in 1:length(S)){
      file <- paste('freqs_M', M[m], '_S', S[s], '.txt', sep = '')
      folder <- './tree2freqs_output/'
      filename <- paste(folder, file, sep = '')
      lines <- read_lines(filename)
      if (length(lines) < 4){next}
      n <- strsplit(lines[3], split = ' #')[[1]][1]
      n <- as.numeric(n)
      n_snv[s,m] <- n
      points(x = M[m], y = n)
    }
  }
  
  num_snv <- c(n_snv[,1], n_snv[,2], n_snv[,3])
  print(length(num_snv))
  m_rate <- c(rep(0.1, 31), rep(1, 31), rep(10, 31))
  snv_df <- data.frame(rate = as.character(m_rate), n_snv = num_snv)
  
  return(snv_df)
}

# MAIN OPERATIONS

rewrite_tree2freq_for_downsample_wrapper()
rewrite_downsample_for_spruce_wrapper()
# snv <- num_snv_analysis()
# ggplot(snv, aes(x=rate, y=log10(n_snv))) + geom_boxplot() + ggtitle('Mutation Count per Rate')
# rewrite_downsample_for_spruce('downsampled_C50_P1_k2_M0.1_S0.txt')
# rewrite_tree2freq_for_downsample_wrapper(M=c(0.1, 1), S = 1:20)
# output <- rewrite_tree2freq_for_downsample('freqs_M0.1_S0.txt')



# output <- rewrite_freq_for_spruce('../frequency output/freq1.txt')
