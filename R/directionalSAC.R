
directionalSAC <- function(community, gradient) {
  
  if (!inherits(community, "data.frame")){
    if(!inherits(community, "matrix"))
      stop("Object community must be of class data.frame or matrix")
  }
  if (!inherits(gradient, "dist")){
    if(!is.vector(gradient)){
     if(!is.matrix(gradient)){
      if(!is.data.frame(gradient))
      stop("Object gradient must be of class vector, matrix, data.frame or dist")
  }}}
  if (any(community < 0))
    stop("Negative value in community")
  if (any(rowSums(community) < 1e-16))
    stop("Empty plots in community")
  if(is.null(rownames(community))) stop("Object community must have row names")
  if(is.vector(gradient)){
    if(length(gradient)!=nrow(community)) stop("Incorrect definition of object gradient")
    if(!is.null(names(gradient))){
    if(any(!rownames(community)%in%names(gradient))) stop("Names in gradient must be the same as row names in community")
    gradient <- gradient[rownames(community)]
    }
    Rgradient <- rownames(community)[order(gradient)]
    agg <- community[Rgradient, ]
    richness <- specnumber(agg)
    average_alfa <- cumsum(richness)/(1:length(richness))
    average <- specaccum(agg, method="collector")$richness
  }
  else {
    if(!is.matrix(gradient)) gradient <- as.matrix(gradient)
    if(ncol(gradient)==1) stop("A gradient with a single quantitative variable must be coded as a vector rather than as a matrix")
    if(nrow(gradient)!=nrow(community)) stop("Incorrect definition of object gradient")
    if(!is.null(rownames(gradient))){
      if(any(!rownames(community)%in%rownames(gradient))) stop("Row names in gradient must be the same as row names in community")
      gradient <- gradient[rownames(community), ]
    }
    res <- array(NA, c(ncol(gradient), nrow(gradient)))
    for(i in 1:ncol(gradient)){
      nami <- rownames(community)
      res[i, ] <- nami[order(gradient[, i])]
    }
    spatial_order <- t(res)
    f <- nrow(spatial_order)
    n <- ncol(spatial_order)
    result <- array(dim = c(f, n))
    alfa_average <- array(dim = c(f, n))
    for(i in 1:n) {
      agg <- community[spatial_order[, i], ]
      richness <- specnumber(agg)
      alfa_s <- cumsum(richness)/(1:length(richness))
      c <- specaccum(agg, method="collector")
      result[, i] <- c$richness
      alfa_average[, i] <- alfa_s
    }
    average <- rowMeans(result)
    average_alfa <- rowMeans(alfa_average)
  }
  beta_s <- average/average_alfa
  beta_S <- (beta_s-1)/((1:length(beta_s))-1)
  exact <- specaccum(community, method = "exact")
  beta_exact <- exact$richness/exact$richness[1]
  beta_N <- (beta_exact-1)/((1:length(beta_exact))-1)
  beta_norm_autocor <- (beta_N- beta_S)/(beta_N+ beta_S)
  SCR <- data.frame(as.matrix(average), as.matrix(exact$richness), as.matrix(average_alfa), as.matrix(beta_s), as.matrix(beta_S), as.matrix(beta_exact), as.matrix(beta_N), as.matrix(beta_norm_autocor))
  names(SCR) <- c("N_SCR", "N_Exact", "Alpha_dir", "Beta_M_dir", "Beta_N_dir","Beta_M", "Beta_N", "Beta_Autocor")
  return(SCR)
}

