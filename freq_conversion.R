# read in frequency output and convert to .data for spruce to enumerate
#!/usr/bin/Rscript

args <- commandArgs(trailingOnly = TRUE)

input_file <- args[1]
output_file <- args[2]
numchain <- args[3]

if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
} else if (length(args)==1) {
  # default output file
  args[2] = "out.txt"
  args[3] = 15
}

library(readr)
library(Canopy)
# source('update_sampling.R')

## CANOPY FUNCTION UPDATE

canopy.sample.nocna = function(R, X, K, numchain, 
                               max.simrun, min.simrun, writeskip, 
                               projectname, cell.line = NULL,
                               plot.likelihood = NULL) {
  if (!is.matrix(R)) {
    stop("R should be a matrix!")
  }
  if (!is.matrix(X)) {
    stop("X should be a matrix!")
  }
  if (min(K) < 2) {
    stop("Smallest number of subclones should be >= 2!\n")
  }
  if (is.null(cell.line)) {
    cell.line = FALSE
  }
  if (is.null(plot.likelihood)) {
    plot.likelihood = TRUE
  }
  if ( plot.likelihood){
    pdf(file = paste(projectname, "_likelihood.pdf", sep = ""), width = 10, height = 5)
  }
  sampname = colnames(R)
  sna.name = rownames(R)
  sampchain = vector("list", length(K))
  ki = 1
  for (k in K) {
    cat("Sample in tree space with", k, "subclones\n")
    sampchaink = vector("list", numchain)
    sampchaink.lik=vector('list',numchain)
    sampchaink.accept.rate=vector('list',numchain)
    for (numi in 1:numchain) {  # numi: number of chain
      cat("\tRunning chain", numi, "out of", numchain, "...\n")
      ###################################### Tree initialization #####
      text = paste(paste(paste(paste("(", 1:(k - 1), ",", sep = ""), 
                               collapse = ""), k, sep = ""), paste(rep(")", (k - 1)), 
                                                                   collapse = ""), ";", sep = "")
      runif.temp=runif(1)
      if(k == 5 & runif.temp<0.5){
        text = c('(1,((2,3),(4,5)));')
      }else if(k == 6 & runif.temp < 1/3){
        text = c('(1,((2,3),(4,(5,6))));')
      }else if(k == 6 & runif.temp > 2/3){
        text = c('(1,(2,((3,4),(5,6))));')
      }else if(k == 7 & runif.temp > 1/4 & runif.temp <= 2/4){
        text=c('(1,((2,3),(4,(5,(6,7)))));')
      }else if(k == 7 & runif.temp > 2/4 & runif.temp <= 3/4){
        text = c('(1,((2,3),((4,5),(6,7))));')
      }else if(k == 7 & runif.temp > 3/4){
        text = c('(1,((2,(3,4)),(5,(6,7))));')
      }
      tree <- read.tree(text = text)
      # print(tree$edge)
      tree$sna = initialsna(tree, sna.name)
      # if(k>=5){tree$relation=getrelation(tree)}
      tree$Z = getZ(tree, sna.name)
      tree$P = initialP(tree, sampname, cell.line)
      tree$VAF = tree$Z%*%tree$P/2
      tree$likelihood = getlikelihood.sna(tree, R, X)
      ###################################### Sample in tree space #####
      sampi = 1
      writei = 1
      samptree = vector("list", max.simrun)
      samptree.lik=rep(NA, max.simrun)
      samptree.accept=rep(NA, max.simrun)
      samptree.accept.rate=rep(NA, max.simrun)
      while(sampi <= max.simrun){
        ######### sample sna positions
        tree.new = tree
        # print(tree$edge)
        tree.new$sna = sampsna(tree)
        # print(tree.new$sna)
        if (any(duplicated(tree.new$sna[,3]))){next}
        tree.new$Z = getZ(tree.new, sna.name)
        # print(tree.new$Z)
        tree.new$VAF = tree.new$Z%*%tree.new$P/2
        tree.new$likelihood = getlikelihood.sna(tree.new, R, X)
        tree.temp=addsamptree(tree,tree.new)
        tree=tree.temp[[1]]
        samptree.accept[sampi]=tree.temp[[2]]
        if (sampi%%writeskip == 0) {
          samptree[[writei]] = tree
          writei = writei + 1
        }
        samptree.lik[sampi]=tree$likelihood
        samptree.accept.rate[sampi]=mean(samptree.accept[max(1,sampi-999):sampi])
        # if ((sampi >= 2*min.simrun) & (samptree.lik[sampi] <= mean(samptree.lik[max((sampi-1000),1):max((sampi-1),1)])) &
        #     (samptree.accept.rate[sampi] <= mean(samptree.accept.rate[max((sampi-1000),1):max((sampi-1),1)])) &
        #     (samptree.accept.rate[sampi] <= 0.1)) break
        sampi = sampi + 1
        ######## sample P (clonal proportions)
        tree.new = tree
        tree.new$P = sampP(tree.new, cell.line)
        tree.new$VAF = tree.new$Z%*%tree.new$P/2
        tree.new$likelihood = getlikelihood.sna(tree.new, R, X)
        tree.temp=addsamptree(tree,tree.new)
        tree=tree.temp[[1]]
        samptree.accept[sampi]=tree.temp[[2]]
        if (sampi%%writeskip == 0) {
          samptree[[writei]] = tree
          writei = writei + 1
        }
        samptree.lik[sampi]=tree$likelihood
        samptree.accept.rate[sampi]=mean(samptree.accept[max(1,sampi-999):sampi])
        # if ((sampi >= 2*min.simrun) & (samptree.lik[sampi] <= mean(samptree.lik[max((sampi-1000),1):max((sampi-1),1)])) &
        #     (samptree.accept.rate[sampi] <= mean(samptree.accept.rate[max((sampi-1000),1):max((sampi-1),1)])) &
        #     (samptree.accept.rate[sampi] <= 0.1)) break
        sampi = sampi + 1
      }
      sampchaink[[numi]] = samptree[1:(writei - 1)]
      sampchaink.lik[[numi]]=samptree.lik
      sampchaink.accept.rate[[numi]]=samptree.accept.rate
    }
    ###################################### plotting and saving #####
    if (plot.likelihood) {
      par(mfrow=c(1,2))
      xmax=ymin=ymax=rep(NA,numchain)
      for(i in 1:numchain){
        xmax[i]=max(which((!is.na(sampchaink.lik[[i]]))))
        ymin[i]=sampchaink.lik[[i]][1]
        ymax[i]=sampchaink.lik[[i]][xmax[i]]
      }
      
      plot(sampchaink.lik[[1]],xlim=c(1,max(xmax)),ylim=c(min(ymin),max(ymax)),
           xlab='Iteration',ylab='Log-likelihood',type='l',
           main=paste('Post. likelihood:',k,'branches'))
      for(numi in 2:numchain){
        points(sampchaink.lik[[numi]],xlim=c(1,max(xmax)),ylim=c(min(ymin),max(ymax)),col=numi,type='l')
      }
      
      plot(sampchaink.accept.rate[[1]],ylim=c(0,1),xlim=c(1,max(xmax)),
           xlab='Iteration',ylab='Acceptance rate',type='l',
           main=paste('Acceptance rate:',k,'branches'))
      for(numi in 2:numchain){
        points(sampchaink.accept.rate[[numi]],ylim=c(0,1),xlim=c(1,max(xmax)),col=numi,type='l')
      }
      par(mfrow=c(1,1))
    }
    sampchain[[ki]] = sampchaink
    ki = ki + 1
  }
  if(plot.likelihood) {
    dev.off()
  }
  return(sampchain)
} 

##
annotate_descendants <- function(tree, root, parents){
  children <- which(tree[root,] == 1)
  tree[parents, children]  <- 2
  
  parents <- c(parents, root)
  for (child in children){
    tree <- annotate_descendants(tree, child, parents)
  }
  return(tree)
}

write_tree <- function(lines, nodes = c()){
  idx <- grep('edges', lines)
  # print(lines)
  num_edges <- as.numeric(strsplit(lines[idx], split = ' #edges')[[1]][1])
  # extract nodes in tree
  if (length(nodes) == 0){
    for (i in (idx+1):(idx+num_edges)){
      if (grepl('_P',lines[i]) | grepl('GL',lines[i])){next}
      elems <- strsplit(lines[i], split = " ")[[1]]
      nodes <- c(nodes, elems)
    }
    nodes <- unique(nodes)
  }
  # print(length(nodes))
  # build tree as matrix; 0 is no relation, 1 is parent, 2 is ancestor; row is higher on tree, column is lower
  tree <- matrix(data = '0', nrow = length(nodes), ncol = length(nodes))
  rownames(tree) <- nodes
  colnames(tree) <- nodes
  
  for (i in (idx+1):(idx+num_edges)){
    if (grepl('_P',lines[i]) | grepl('GL',lines[i])){next}
    elems <- strsplit(lines[i], split = " ")[[1]]
    tree[elems[1], elems[2]] <- 1
    # if (all(elems %in% rownames(tree))){
    #   tree[elems[1], elems[2]] <- 1
    # }
  }
  
  tree <- annotate_descendants(tree, 1, c())
  return(tree)
}

compare_trees <- function(target_tree, ref_tree){
  
  branching_mtx <- function(tree){
    mtx <- tree != 0
    mtx[lower.tri(mtx, diag = TRUE)] <- TRUE
    # paths <- which(mtx)
    # mtx <- t(mtx)
    # mtx[paths] <- TRUE
    return(mtx)
  }
  
  # parent == 1, ancecstor == 2
  true_parents <- length(which(ref_tree == 1))
  true_ancestors <- length(which(ref_tree == 1 | ref_tree == 2))

  captured_parents <- length(which(target_tree == 1 & ref_tree == 1))
  captured_ancestors <- length(which(target_tree != 0 & ref_tree != 0))
  
  parentage_accuracy <- captured_parents/true_parents
  ancestor_accuracy <- captured_ancestors/true_ancestors
  
  target_branching <- branching_mtx(target_tree)
  ref_branching <- branching_mtx(ref_tree)
  
  true_branching <- length(which(!ref_branching))
  captured_branching <- length(which(!target_branching & !ref_branching))
  
  branching_accuracy <- captured_branching/true_branching
    
  return(c(parentage_accuracy, ancestor_accuracy, branching_accuracy))
}

iterate_inferred_tree_list <- function(inferred_lines, true_lines){
  true_tree <- write_tree(true_lines)
  nodes <- rownames(true_tree)
  idxs <- grep('#edges', inferred_lines)
  n_trees <- length(idxs)
  # print(n_trees)
  accuracy <- matrix(data = 0, nrow = n_trees, ncol = 3)
  colnames(accuracy) <- c('direct path', 'indirect path', 'no path')
  accuracy_idx <- 1
  
  for (i in idxs){
    num_edges <- as.numeric(strsplit(inferred_lines[i], split = ' #edges')[[1]][1])
    # print(num_edges)
    tree <- write_tree(inferred_lines[i:(i+num_edges)], nodes = nodes)
    compar <- compare_trees(tree, true_tree)
    accuracy[accuracy_idx,1] <- compar[1]
    accuracy[accuracy_idx,2] <- compar[2]
    accuracy[accuracy_idx,3] <- compar[3]
    accuracy_idx <- accuracy_idx+1
  }
  
  return(accuracy)
}

compare_spruce_output <- function(inferred_tree_file, true_tree_file, output_filename){
  true_tree <- read_lines(true_tree_file)
  inferred_trees <- read_lines(inferred_tree_file)
  # n_inferred_lines <- length(inferred_trees)
  # inferred_trees <- inferred_trees[2:n_inferred_lines]
  
  accuracy_matrix <- iterate_inferred_tree_list(inferred_trees, true_tree)
  write.table(accuracy_matrix, file = output_filename)
  return(accuracy_matrix)
}

rewrite_tree2freq_for_downsample <- function(filename, coverage = 10000, input_folder = '~/GitHub/VAFFP-NonUniq/tree2freqs_output/', 
                                             output_folder = '~/GitHub/VAFFP-NonUniq/downsample_input/'){
  output_file <- paste(output_folder, filename, sep = '')
  file <- read_lines(paste(input_folder, filename, sep = ''))
  if (length(file) < 2){return()}
  # file <- file[-c(1)]
  
  table <- file[c(5:length(file))]
  write_lines(table, 'temp.txt')
  mtx <- read_delim("temp.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  colnames(mtx) <- c('sample_index','sample_label','anatomical_site_index','anatomical_site_label',
                     'character_index','character_label','ref','var')
  
  mtx[,7] <- floor(coverage*(1-mtx[,8]))
  mtx[,8] <- coverage-mtx[,7]
  # mtx <- format(mtx, scientific = FALSE)
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

rewrite_downsample_for_spruce <- function(freqs_filename, downsample_filename, coverage = 10000, 
                                          downsample_input_folder = '~/GitHub/VAFFP-NonUniq/downsample_output/', 
                                          freqs_input_folder = '~/GitHub/VAFFP-NonUniq/reads2freqs_output/', 
                                          output_folder = '~/GitHub/VAFFP-NonUniq/spruce_input/'){
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
  colnames(downsample_mtx) <- c('sample_index','sample_label','anatomical_site_index',
                                'anatomical_site_label','character_index','character_label','ref','var')
  f_avg <- downsample_mtx[,8]/(downsample_mtx[,8]+downsample_mtx[,7])
  
  table <- freqs_file[c(5:length(freqs_file))]
  write_lines(table, 'temp.txt')
  freqs_mtx <- read_delim("temp.txt", "\t", escape_double = FALSE, col_names = FALSE, trim_ws = TRUE)
  colnames(freqs_mtx) <- c('sample_index','sample_label','anatomical_site_index','anatomical_site_label',
                           'character_index','character_label','f-','f+')
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
  output_mtx[,5] <- f_minus/2
  output_mtx[,6] <- f_avg/2
  output_mtx[,7] <- f_plus/2
  
  # set x,y,mu
  output_mtx[,8] <- as.integer(1)
  output_mtx[,9] <- as.integer(1)
  output_mtx[,10] <- as.integer(1)
  
  colnames(output_mtx) <- c('sample_index','sample_label','character_index','character_label','f-','f','f+','x','y','mu')
  
  write_delim(output_mtx, 'temp.txt', delim = " ")
  new_file <- read_lines('temp.txt')
  new_file[1] <- paste('#', new_file[1], sep = '')
  # print(file[1:2])
  new_file <- c(freqs_file[2:3], new_file)
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
      folder <- '~/GitHub/VAFFP-NonUniq/tree2freqs_output/'
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

iterate_accuracy <- function(S = 0:30, M = c(0.1,0.2,0.3,0.4), K = c(1,2,5,10)){
  input_folder <- 'enumerate_output/'
  output_folder <- 'accuracy_output/'
  base_input_filename <- '_clustered.txt'
  base_output_filename <- '_accuracy.txt'
  
  tree_folder <- 'simulate_output/'
  tree_base_filename <- '.tree'
  
  for (m in M){
    print(m)
    for (s in S){
      print(s)
      tree_filename <- paste(tree_folder,'M',m,'/T_seed',s,tree_base_filename, sep = '')
      for (k in K){
        print(k)
        input_filename <- paste(input_folder,'M',m,'_S',s,'_k',k,base_input_filename, sep = '')
        if (file.exists(input_filename)){
          output_filename <- paste(output_folder,'M',m,'_S',s,'_k',k,base_output_filename, sep = '')
          inferred_output <- compare_spruce_output(input_filename, tree_filename, output_filename)
        }
      }
    }
  }
}

extract_cluster_frequency <- function(filename){
  lines <- read_lines(filename)
  info_lines <- grep('#', lines)
  table <- lines[c((length(info_lines)+1):length(lines))]
  write_lines(table, 'temp.txt')
  mtx <- read_delim("temp.txt", "\t", escape_double = FALSE, 
                    col_names = c('sample_index',	'sample_label',	'anatomical_site_index',	'anatomical_site_label',	
                                  'character_index',	'character_label',	'f-',	'f+'),
                    trim_ws = TRUE)
  return(mtx)
}

construct_canopy_input <- function(mtx, nReads = 10000){

  sample_label <- unlist(unique(mtx[,'sample_label']))
  character_label <- unlist(unique(mtx[,'character_label']))
  
  k <- length(sample_label)
  nSnv <- length(character_label)
  
  var_reads <- unlist(mtx[,'f-'])*nReads
  
  chrom_names <- paste('Chrom_', seq(1,nSnv,1), sep = '')
  nChrom <- nSnv
  nM <- 2 # major chromosomes
  nm <- 0 # minor chromosomes
  epsilon <- 0
  
  WM <- matrix(data = nM, nrow = nChrom, ncol = k, dimnames = list(chrom_names, sample_label))
  Wm <- matrix(data = nm, nrow = nChrom, ncol = k, dimnames = list(chrom_names, sample_label))

  # Y <- matrix(data = 1, nrow = nSnv, ncol = nChrom+1, dimnames = list(character_label, c('non-cna_region', chrom_names)))
  Y <- diag(x = 1, nrow = nSnv, ncol = nSnv)
  Y <- cbind(0, Y)
  rownames(Y) <- character_label
  colnames(Y) <- c('non-cna_region', chrom_names)
  # Y[,1] <- 0
  
  R <- matrix(data = var_reads, nrow = nSnv, ncol = k, dimnames = list(character_label, sample_label))
  X <- matrix(data = nReads, nrow = nSnv, ncol = k, dimnames = list(character_label, sample_label))
  # rownames(R) <- character_label
  # colnames(R) <- sample_label
  
  list(WM = WM, Wm = Wm, Y = Y, R = R, X = X, epsilon = epsilon, K = nSnv)
  # list(R = R, X = X, K = nSnv)
}

write_canopy_trees <- function(trees_list){
  get_edges <- function(tree, tree_idx, lines){
    edges <- tree$edge
    new_lines <- c()
    edge <- tree$sna[,3]
    muts <- as.vector(rownames(tree$sna))
    map <- cbind(muts, edge)
    for (i in 1:nrow(edges)){
      if (!(edges[i,1] %in% map[,2])){next}
      if (!(edges[i,2] %in% map[,2])){
        new_pair <- c(map[which(map[,2] == edges[i,1]),1], edges[i,2])
        map <- rbind(map, new_pair)
      }
      idx1 <- which(map[,2] == edges[i,1])
      idx2 <- which(map[,2] == edges[i,2])
      mut1 <- map[idx1,1]
      mut2 <- map[idx2,1]
      if (mut1 == mut2){next}
      new_line <- paste(mut1, mut2, sep = ' ')
      new_lines <- c(new_lines, new_line)
    }
    
    new_lines <- c(paste(length(new_lines), '#edges, tree', tree_idx-1, sep = ' '), new_lines)
    lines <- c(lines, new_lines)
    return(lines)
  }
  
  num_trees <- length(trees_list)
  
  lines <- c(paste(num_trees, '#trees', sep = ' '))
  for (i in 1:length(trees_list)){
    lines <- get_edges(trees_list[[i]], i, lines)
  }
  
  return(lines)
}

# MAIN OPERATIONS

# rewrite_tree2freq_for_downsample_wrapper()
# rewrite_downsample_for_spruce_wrapper()
# snv <- num_snv_analysis()
# ggplot(snv, aes(x=rate, y=log10(n_snv))) + geom_boxplot() + ggtitle('Mutation Count per Rate')
# rewrite_downsample_for_spruce('downsampled_C50_P1_k2_M0.1_S0.txt')
# rewrite_tree2freq_for_downsample_wrapper(M=c(0.1, 1), S = 1:20)
# output <- rewrite_tree2freq_for_downsample('freqs_M0.1_S0.txt')

# lines <- read_lines('~/GitHub/VAFFP-NonUniq/simulate_output/M0.1/T_seed0.tree')
# tree <- write_tree(lines)
# 
# inferred_output <- compare_spruce_output('~/GitHub/VAFFP-NonUniq/enumerate_output/M0.1_S0_k1_clustered.txt', 
#                                               '~/GitHub/VAFFP-NonUniq/simulate_output/M0.1/T_seed0.tree', 'test.txt')
# iterate_accuracy()
# output <- rewrite_freq_for_spruce('../frequency output/freq1.txt')
projectname <- 'canopy_sampling'
# input_file <- paste('mix_output/', projectname, '_clustered.tsv', sep = '')

mtx <-extract_cluster_frequency(input_file)
canopy_input <- construct_canopy_input(mtx, nReads = 10^10)

# test_canopy_output <- canopy.sample(R = canopy_input$R, X = canopy_input$X, WM = canopy_input$WM, Wm = canopy_input$Wm, 
#                                     epsilonM = canopy_input$epsilon, epsilonm = canopy_input$epsilon, Y = canopy_input$Y)
R <- canopy_input$R
X <- canopy_input$X
K <- canopy_input$K+1
# numchain <- 2
# tic()
sampchain <- canopy.sample.nocna(R = R, X = X, K = K,
                                 # WM = canopy_input$WM, Wm = canopy_input$Wm, epsilonM = 0.001, epsilonm = 0.001,Y = canopy_input$Y,
                                 numchain = numchain, max.simrun = 100000,
                                 min.simrun = 20000, writeskip = 200,
                                 projectname = projectname, cell.line = TRUE, plot.likelihood = TRUE)

burnin <- 100
thin <- 5
post <- canopy.post(sampchain = sampchain, projectname = projectname, K = K,
                   numchain = numchain, burnin = burnin, thin = thin, 
                   optK = K, post.config.cutoff = 0.05)
samptreethin = post[[1]]   # list of all post-burnin and thinning trees

# no_clustering <- c()
# for (i in 1:length(sampchain[[1]][[1]])){
#   if (!any(duplicated(sampchain[[1]][[1]][[i]]$sna[,3]))){
#     no_clustering <- c(no_clustering, i)
#   }
# }
# # toc()
# # for (i in 1:length(sampchain[[1]][[2]])){
# #   if (!any(duplicated(sampchain[[1]][[2]][[i]]$sna[,3]))){
# #     no_clustering <- c(no_clustering, i)
# #   }
# # }
# print(length(no_clustering))
# print(length(sampchain[[1]][[1]]))
# no_clustering_trees <- sampchain[[1]][[1]][no_clustering]
lines <- write_canopy_trees(samptreethin)
write_lines(lines, output_file)
# print('finished')
# canopy.plottree(samptreethin[[1]], pdf = TRUE, pdf.name = 'tree.pdf')
# for (i in 1:5){
#   filename <- paste('tree_', i, '.pdf', sep = '')
#   canopy.plottree(canopy_output[[1]][[2]][[i]], pdf = TRUE, pdf.name = filename)
# }
