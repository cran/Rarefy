
rao_permuted<-function(comm_st,dist_f,random=99) {
  
  if (!inherits(comm_st, "matrix") && !inherits(comm_st, "data.frame")) stop("Non convenient comm_st")
  if (any(comm_st<0)) stop("Negative value in comm_st")
  if(suppressWarnings(any(rowSums(comm_st))==0)) stop("Empty row")
  if(suppressWarnings(any(colSums(comm_st))==0)) {
    v<-apply(comm_st, 2, function(col) any(col !=0 ))
    comm_st<-comm_st[,v]
  }
  
  if (!inherits(dist_f, "dist")) stop("Object of class 'dist' expected for dist_f")
  if (!is.euclid(sqrt(dist_f))) warning("Squared Euclidean or Euclidean property expected for dist_f")
  dist_f<-as.matrix(dist_f)
  if (ncol(comm_st) >= nrow(dist_f)) stop("comm_st must have a minor number of species than dist_f")
  if (any(dist_f<0)) stop("Negative value in dist_f")
  
  comm_st<-sweep(comm_st,1,rowSums(comm_st),"/")
  nami<-colnames(comm_st)
  r_fin<-array(dim = c(random,nrow(comm_st)))
  for(i in 1:random) {
    r<-sample(1:nrow(comm_st),nrow(comm_st))
    s<-sample(1:nrow(dist_f),ncol(comm_st))
    x<-comm_st[r,]
    disfx<-dist_f[s,s]
    x<-apply(x,2,cummean)
    r_fin[i,]<-apply(x,1,function(newp) t(as.vector(newp)) %*% disfx %*% as.vector(newp))
  }
  
  rare<-colMeans(r_fin)
  IC_up <- rare + (1.96*(sd(r_fin)/sqrt(random)))
  IC_low <- rare - (1.96*(sd(r_fin)/sqrt(random)))
  df<-data.frame(rare,IC_up,IC_low)
  colnames(df)<-c('Rao','IC_up','IC_low')
  
  return(df)
  
}

