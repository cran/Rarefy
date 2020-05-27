
rarexpected_fun<-function(comm,dist_f) {
  
  if(is.null(colnames(comm))) stop("comm must have names for columns") 
  if (!inherits(comm, "matrix") && !inherits(comm, "data.frame")) stop("Non convenient comm")
  if (any(comm<0)) stop("Negative value in comm")
  if(suppressWarnings(any(rowSums(comm))==0)) stop("Empty row")
  if(!all(apply(comm,1, function(x) all(x %in% c(0,1))))) stop("Non convenient comm")
  if(suppressWarnings(any(colSums(comm))==0)) {
    v<-apply(comm, 2, function(col) any(col !=0 ))
    comm<-comm[,v]
  }
  
  if (!inherits(dist_f, "dist")) stop("Object of class 'dist' expected for dist_f")
  if (!is.euclid(sqrt(dist_f))) warning("Squared Euclidean or Euclidean property expected for dist_f")
  dist_f<-as.matrix(dist_f)
  if (ncol(comm) != nrow(dist_f)) stop("comm and dist_f don't have the same number of species")
  if(!is.null(rownames(dist_f)) && !is.null(colnames(comm)) && any(!colnames(comm)%in%rownames(dist_f))) stop("Names in dist_f must be the same as col names in comm")
  if (any(dist_f<0)) stop("Negative value in dist_f")
  
  freq <- colSums(comm > 0)
  fin<-rep(NA,nrow(comm))
  for (i in 1:nrow(comm)) {
    v<-1-ifelse(nrow(comm) - freq < i, 0 , exp(lchoose(nrow(comm) - freq,i)-lchoose(nrow(comm),i)))
    v<-as.matrix(v/sum(v))
    if (sum(v) < 1e-16) 
      fin[i] <- 0
    else fin[i] <- (t(v) %*% (as.matrix(dist_f)) %*% v) 
  }
  
  return(fin)
  
}


