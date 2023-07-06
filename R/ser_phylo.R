
ser_phylo <- function(comm, tree = NULL, dist_xy, method = c("faith", "barker", "Ia", "hill", "tsallis", "renyi", "fun_div"), exp = 0, fun_div = NULL, args = NULL, verbose = FALSE, comparison = FALSE, resampling = 99) {
  
  method <- method[1]
  if(!method %in% c("faith","barker","Ia","tsallis", "hill", "renyi", "fun_div")) stop("Unavailable method")
  
  if(is.null(colnames(comm))) stop("comm must have names for columns") 
  if(!inherits(comm, "matrix") && !inherits(comm, "data.frame")) stop("Non convenient comm")
  if(any(comm < 0)) stop("Negative value in comm")
  if(suppressWarnings(any(rowSums(comm)) == 0)) stop("Empty row")
  if(suppressWarnings(any(colSums(comm)) == 0)) {
    v <- apply(comm, 2, function(col) any(col != 0))
    comm <- comm[,v]
  }
  
  if (!inherits(dist_xy, "dist")){
    if(!is.vector(dist_xy)){
      if(!is.matrix(dist_xy)){
        if(!is.data.frame(dist_xy))
          stop("Object dist_xy must be of class vector, matrix, data.frame or dist") }}}
  if(!inherits(dist_xy, "dist")) {
    if(is.vector(dist_xy)){
      if(length(dist_xy) != nrow(comm)) stop("Incorrect definition of object gradient")
      if(!is.null(names(dist_xy))){
        if(any(!rownames(comm) %in% names(dist_xy))) stop("Names in gradient must be the same as row names in community")
        dist_xy <- dist_xy[rownames(comm)] }}
    else { 
      if (nrow(comm) != nrow(dist_xy)) stop("comm and dist_xy don't have the same number of plots") 
      }
    dist_xy <- dist(scale(dist_xy)) 
    warning("dist_xy is used to calculate a distance matrix between sampling units using the method euclidean of the function dist()")
  }
  dist_xy <- as.matrix(dist_xy)
  if (nrow(comm) != nrow(dist_xy)) stop("comm and dist_xy don't have the same number of plots")
  if (any(dist_xy<0)) stop("Negative value in dist_xy")
  if(!is.null(rownames(dist_xy)) && !is.null(rownames(comm)) ) {
    if(any(!rownames(comm) %in% rownames(dist_xy))) stop("comm and dist_xy must have the same names for the plots")
  } else if(!is.null(rownames(dist_xy)) && is.null(rownames(comm))) {
    rownames(comm) <- rownames(dist_xy)
    warning("comm has no row names")
    warning("row names of dist_xy set as row names of comm")
  } 
  
  if(is.null(tree) && method != 'fun_div') stop("tree must have a value")
  if(!is.null(tree)) {
    if(!inherits(tree, "phylo") && !inherits(tree, "phylo4")) stop("tree must be of class phylo or phylo4") 
    if(inherits(tree, "phylo4")) suppressWarnings(tree <- as(tree, "phylo"))
    if(any(!colnames(comm) %in% tree$tip.label)) stop("comm contains tip names that are not available in tree")
    if(is.null(tree$edge.length)) stop("edge lengths are required for tree") 
  }
  
  if(method == 'faith') {
    r_fin <- array(dim = c(ncol(dist_xy), ncol(dist_xy)))
    nami <- rownames(comm)
    for(i in 1:ncol(dist_xy)) {
      v <- nami[order(dist_xy[, i])]
      x <- comm[v,]
      x <- apply(x, 2, cumsum)
      sub_tree <- apply(x, 1, function(x) drop.tip(tree,which(!(tree$tip.label %in% colnames(comm)[x > 0]))))
      r_fin[i,] <- unlist(lapply(sub_tree,function(x) sum(x$edge.length)))
    }
    rare<-colMeans(r_fin)
  }
  
  else if(method == 'barker') {
    comm <- sweep(comm, 2, colSums(comm), "/")
    r_fin <- array(dim = c(ncol(dist_xy), ncol(dist_xy)))
    nami <- rownames(comm)
    
    get.leaves <- function(x, st){
      leaves.node <- tips(st, x[2])
    }
    
    for(i in 1:ncol(dist_xy)) {
      v <- nami[order(dist_xy[, i])]
      x <- comm[v,]
      x <- apply(x, 2, cumsum)
      sub_tree <- apply(x, 1, function(x) drop.tip(tree,which(!(tree$tip.label %in% colnames(comm)[x > 0]))))
      lv <- lapply(sub_tree, function(y) apply(y$edge, 1, function(z) get.leaves(z,y)))
      for(j in 1:nrow(x)) {
        v <- x[j,]
        rel <- unlist(lapply(lv[[j]], function(l) mean(v[l])))
        r_fin[i,j] <- nrow(sub_tree[[j]]$edge) * ((sum(sub_tree[[j]]$edge.length * rel)) / sum(rel))
      }
    }
    rare <- colMeans(r_fin)
  }
  
  else if(method == 'Ia') {
    r_fin <- array(dim = c(ncol(dist_xy), ncol(dist_xy)))
    nami <- rownames(comm)
    for(i in 1:ncol(dist_xy)) {
      v <- nami[order(dist_xy[, i])]
      x <- comm[v,]
      x <- apply(x, 2, cumsum)
      r_fin[i,] <- pIa(tree, x, exponent = exp)[[1]]
    }
    rare <- colMeans(r_fin)
  }
  
  else if(method == 'hill' || method == 'tsallis' || method == 'renyi') {
    r_fin <- array(dim = c(ncol(dist_xy), ncol(dist_xy)))
    nami <- rownames(comm)
    for(i in 1:ncol(dist_xy)) {
      v <- nami[order(dist_xy[, i])]
      x <- comm[v,]
      x <- apply(x, 2, cumsum)
      r_fin[i,] <- evodivparam(tree, x, method = method, q = exp)
    }
    rare <- colMeans(r_fin)
  }
  
  else if(method == 'fun_div') {
    if(!inherits(fun_div, 'character')) stop("fun_div must be a character") 
    if(!exists(fun_div)) stop("the function doesn't exist") 
    
    f <- match.fun(fun_div)
    arg <- as.list(args(f)) 
    v <- names(arg)
    
    if(verbose) {
      v[length(v)] <- 'stop'
      ch <- 1
      i <- 2
      l <- list(comm)
      cat('Wich argument is the community matrix?')
      cat(paste(1:(length(v)-1), '-', v[-length(v)]))
      ch <- readline("Argument number: ")
      ch <- as.numeric(ch)
      if(ch %in% 1:length(v) && !v[ch] %in% c('stop')) n <- v[ch] else stop(paste(ch, "is not a possible choice"))
      cat('Which arguments do you want to set?')
      
      while(!v[ch] %in% c('stop')) {
        cat(paste(1:length(v), '-', v)) 
        ch<-readline("Argument number: ") 
        ch<-as.numeric(ch)
        if(v[ch] == 'stop') break()
        if(ch %in% 1:length(v) && !v[ch] %in% c('stop')) {
          n[i] <- v[ch]
          v[ch] <- readline(paste(v[ch], ' = '))
        }
        if(exists(v[ch]) && v[ch] != n[i]) { 
          l[[i]] <- get(v[ch])
          i <- i+1
        }
        else if(v[ch] %in% c('FALSE', 'TRUE', 'T', 'F') && v[ch] != n[i]) {
          l[[i]] <- as.logical(v[ch])
          i <- i+1
        }
        else if(!is.na(as.numeric(v[ch])) && v[ch]!= n[i]) {
          l[[i]]<-as.numeric(v[ch])
          i <- i+1
        }
        else if(grepl('^c\\(.+\\)$', v[ch]) && v[ch] != n[i]) {
          v[ch] <- gsub('c\\(', '', v[ch])
          v[ch] <- gsub('\\)', '', v[ch])
          v[ch] <- gsub('\'', '', v[ch])
          st <- strsplit(v[ch], ',')
          if(all(suppressWarnings(!is.na(sapply(st, as.numeric))))) {
            st <- unlist(lapply(st, as.numeric))
            l[[i]] <- st
            i <- i+1
          }
          else if(all(st %in% c('FALSE', 'TRUE', 'T', 'F'))) {
            st <- unlist(lapply(st, as.logical))
            l[[i]] <- st
            i <- i+1
          }
          else {
            l[[i]] <- st
            i <- i+1
          }
        }
        else if(v[ch] != n[i]) {
          l[[i]] <- v[ch]
          i <- i+1
        }
      }
      
      names(l) <- n
      
      nami <- rownames(comm)
      r_fin <- array(dim = c(nrow(comm), nrow(comm)))
      for(i in 1:nrow(comm)) {
        v <- nami[order(dist_xy[, i])]
        x <- comm[v,]
        x <- apply(x,2,cumsum)
        l[[1]] <- x
        r_fin[i,] <- do.call(f, l)
      }
      rare<-colMeans(r_fin)
    }
    
    else {
      if(!all(names(args) %in% v)) stop("The arguments must be the ones specified by the function")
      
      ind <- match(NA, unlist(args))
      ind <- match(names(unlist(args)[ind]), names(args))
      
      nami <- rownames(comm)
      r_fin <- array(dim = c(nrow(comm), nrow(comm)))
      for(i in 1:nrow(comm)) {
        v <- nami[order(dist_xy[, i])]
        x <- comm[v,]
        x <- apply(x, 2, cumsum)
        args[[ind]] <- x
        r_fin[i,] <- do.call(f, args)
      }
      rare <- colMeans(r_fin)
    }
  }
  
  if(comparison) {
    rare1 <- rare_phylo(comm, tree, method = method, exp = exp, resampling = resampling, fun_div = fun_div, args = args)
    IC_up <- rare + (1.96 * (sd(r_fin) / sqrt(nrow(comm))))
    IC_low <- rare - (1.96 * (sd(r_fin) / sqrt(nrow(comm))))
    df <- data.frame(rare, IC_up, IC_low, rare1)
    colnames(df) <- c('Ser', 'IC_up_spaz', 'IC_low_spaz', 'Rarefaction', 'IC_up', 'IC_low')
  }
  
  else {
    IC_up <- rare + (1.96 * (sd(r_fin) / sqrt(nrow(comm))))
    IC_low <- rare - (1.96 * (sd(r_fin) / sqrt(nrow(comm))))
    df <- data.frame(rare, IC_up, IC_low)
    colnames(df) <- c('Rarefaction', 'IC_up', 'IC_low')
  }
  
  return(df)
  
}

