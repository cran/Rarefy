
rare_beta<-function(comm,dist_xy=NULL,method=c("whittaker","jaccard","bray","cody","fun_div"),random=99,fun_div=NULL,args=NULL,verbose=FALSE,spatial=FALSE) {
  
  method <- method[1]
  if(!method%in%c("whittaker","jaccard","sorensen","bray","cody","fun_div")) stop("Unavailable method")
  if(method=='cody' && !spatial) stop("cody's index can be used only with spatial=TRUE")
  
  if(is.null(colnames(comm))) stop("comm must have names for columns") 
  if (!inherits(comm, "matrix") && !inherits(comm, "data.frame")) stop("Non convenient comm")
  if (any(comm<0)) stop("Negative value in comm")
  if(suppressWarnings(any(rowSums(comm))==0)) stop("Empty row")
  if(suppressWarnings(any(colSums(comm))==0)) {
    v<-apply(comm, 2, function(col) any(col !=0 ))
    comm<-comm[,v]
  }
  
  
  if(!is.null(dist_xy)) {
    if (!inherits(dist_xy, "dist")) stop("Object of class 'dist' expected for dist_xy") 
    dist_xy<-as.matrix(dist_xy)
    if (nrow(comm) != nrow(dist_xy)) stop("comm and dist_xy don't have the same number of plots")
    if(!is.null(rownames(dist_xy)) && !is.null(rownames(comm)) ) {
      if(any(!rownames(comm)%in%rownames(dist_xy))) stop("comm and dist_xy must have the same names for the plots")
    } else if(!is.null(rownames(dist_xy)) && is.null(rownames(comm))) {
      rownames(comm)<-rownames(dist_xy)
      warning("comm has no row names")
      warning("row names of dist_xy set as row names of comm")
    } 
    if (any(dist_xy<0)) stop("Negative value in dist_xy") 
  }
  if(spatial && is.null(dist_xy)) stop("dist_xy is requiered for the spatial explicit rarefaction")
  
  if(method=='whittaker') {
    r_fin<-array(dim = c(ifelse(spatial,nrow(comm),random),nrow(comm)))
    nami<-rownames(comm)
    for(i in 1:ifelse(spatial,nrow(comm),random)) {
      if(spatial) v<-nami[order(dist_xy[,i])]
      else v<-sample(1:nrow(comm),nrow(comm))
      x<-comm[v,]
      sub<-specnumber(x)
      sub<-cummean(sub)
      x<-apply(x,2,cumsum)
      g<-specnumber(x)
      r_fin[i,]<-(g-sub)/sub
    }
    rare<-colMeans(r_fin)
    IC_up <- rare + (1.96*(sd(r_fin)/sqrt(ifelse(spatial,nrow(comm),random))))
    IC_low <- rare - (1.96*(sd(r_fin)/sqrt(ifelse(spatial,nrow(comm),random))))
  }
  
  else if(method=='jaccard' || method=='bray') {
    r_fin<-array(dim = c(ifelse(spatial,nrow(comm),random),nrow(comm)-1))
    nami<-rownames(comm)
    for(i in 1:ifelse(spatial,nrow(comm),random)) {
      if(spatial) v<-nami[order(dist_xy[,i])]
      else v<-sample(1:nrow(comm),nrow(comm))
      x<-comm[v,]
      for(j in 2:nrow(comm)) {
        sub<-x[1:j,]
        r_fin[i,(j-1)]<-mean(as.matrix(vegdist(sub,method = method)))
      }
    }
    rare<-colMeans(r_fin)
    IC_up <- rare + (1.96*(sd(r_fin)/sqrt(ifelse(spatial,nrow(comm),random))))
    IC_low <- rare - (1.96*(sd(r_fin)/sqrt(ifelse(spatial,nrow(comm),random))))
    rare<-c(NA,rare)
    IC_up<-c(NA,IC_up)
    IC_low<-c(NA,IC_low)
  }
  
  
  else if(method=='cody') {
    r_fin<-array(dim = c(nrow(comm),nrow(comm)-1))
    nami<-rownames(comm)
    for(i in 1:nrow(comm)) {
      v<-nami[order(dist_xy[,i])]
      x<-comm[v,]
      p1<-x[1,]
      x1<-apply(x,2,cumsum)
      g<-apply(x1[-1,],1,function(x) length(x[x[p1==0]>0]))
      l<-apply(x[-1,],1,function(x) length(x[x==0]))
      r_fin[i,]<-(g+l)/2
    }
    rare<-colMeans(r_fin)
    IC_up <- rare + (1.96*(sd(r_fin)/sqrt(ifelse(spatial,nrow(comm),random))))
    IC_low <- rare - (1.96*(sd(r_fin)/sqrt(ifelse(spatial,nrow(comm),random))))
    rare<-c(NA,rare)
    IC_up<-c(NA,IC_up)
    IC_low<-c(NA,IC_low)
  }
  
  else if(method=='fun_div') {
    if(!inherits(fun_div, 'character')) stop("fun_div must be a character") 
    if(!exists(fun_div)) stop("the function doesn't exist") 
    
    f<-match.fun(fun_div)
    arg<-as.list(args(f)) 
    v<-names(arg)
    
    if(verbose) {
      v[length(v)]<-'stop'
      ch<-1
      i<-2
      l<-list(comm)
      cat('Wich argument is the community matrix?')
      cat(paste(1:(length(v)-1),'-',v[-length(v)]))
      ch<-readline("Argument number: ")
      ch<-as.numeric(ch)
      if(ch %in% 1:length(v) && !v[ch] %in% c('stop')) n<-v[ch] else stop(paste(ch,"is not a possible choice"))
      cat('Which arguments do you want to set?')
      
      while(!v[ch] %in% c('stop')) {
        cat(paste(1:length(v),'-',v)) 
        ch<-readline("Argument number: ") 
        ch<-as.numeric(ch)
        if(v[ch]=='stop') break()
        if(ch %in% 1:length(v) && !v[ch] %in% c('stop')) {
          n[i]<-v[ch]
          v[ch]<-readline(paste(v[ch], ' = '))
        }
        if(exists(v[ch]) && v[ch]!= n) { 
          l[[i]]<-get(v[ch])
          i<-i+1
        }
        else if(v[ch] %in% c('FALSE','TRUE','T','F') && v[ch]!= n) {
          l[[i]]<-as.logical(v[ch])
          i<-i+1
        }
        else if(!is.na(as.numeric(v[ch])) && v[ch]!= n) {
          l[[i]]<-as.numeric(v[ch])
          i<-i+1
        }
        else if(grepl('^c\\(.+\\)$',v[ch]) && v[ch]!= n[i]) {
          v[ch]<-gsub('c\\(','',v[ch])
          v[ch]<-gsub('\\)','',v[ch])
          v[ch]<-gsub('\'','',v[ch])
          st<-strsplit(v[ch],',')
          if(all(suppressWarnings(!is.na(sapply(st,as.numeric))))) {
            st<-unlist(lapply(st,as.numeric))
            l[[i]]<-st
            i<-i+1
          }
          else if(all(st %in% c('FALSE','TRUE','T','F'))) {
            st<-unlist(lapply(st,as.logical))
            l[[i]]<-st
            i<-i+1
          }
          else {
            l[[i]]<-st
            i<-i+1
          }
        }
        else if(v[ch]!= n[i]) {
          l[[i]]<-gsub("\'|\"",'',v[ch])
          i<-i+1
        }
      }
      
      names(l)<-n
      
      r_fin<-array(dim = c(ifelse(spatial,nrow(comm),random),nrow(comm)-1))
      nami<-rownames(comm)
      for(i in 1:ifelse(spatial,nrow(comm),random)) {
        if(spatial) v<-nami[order(dist_xy[,i])]
        else v<-sample(1:nrow(comm),nrow(comm))
        x<-comm[v,]
        for(i in 2:nrow(comm)) {
          l[[1]]<-x[1:j,]
          r_fin[i,]<-do.call(f,l)
        }
      }
      rare<-colMeans(r_fin)
      IC_up <- rare + (1.96*(sd(r_fin)/sqrt(ifelse(spatial,nrow(comm),random))))
      IC_low <- rare - (1.96*(sd(r_fin)/sqrt(ifelse(spatial,nrow(comm),random))))
      rare<-c(NA,rare)
      IC_up<-c(NA,IC_up)
      IC_low<-c(NA,IC_low)
    }
    
    else {
      if(!all(names(args) %in% v)) stop("The arguments must be the ones specified by the function")
      
      ind<-match(NA,unlist(args))
      ind<-match(names(unlist(args)[ind]), names(args))
      
      r_fin<-array(dim = c(ifelse(spatial,nrow(comm),random),nrow(comm)-1))
      nami<-rownames(comm)
      for(i in 1:ifelse(spatial,nrow(comm),random)) {
        if(spatial) v<-nami[order(dist_xy[,i])]
        else v<-sample(1:nrow(comm),nrow(comm))
        x<-comm[v,]
        for(i in 2:nrow(comm)) {
          l[[1]]<-x[1:j,]
          r_fin[i,]<-do.call(f,l)
        }
      }
    }
    rare<-colMeans(r_fin)
    IC_up <- rare + (1.96*(sd(r_fin)/sqrt(ifelse(spatial,nrow(comm),random))))
    IC_low <- rare - (1.96*(sd(r_fin)/sqrt(ifelse(spatial,nrow(comm),random))))
    rare<-c(NA,rare)
    IC_up<-c(NA,IC_up)
    IC_low<-c(NA,IC_low)
  }
  
  df<-data.frame(rare,IC_up,IC_low)
  colnames(df)<-c('Rarefaction','IC_up','IC_low')
  
  return(df)
  
}

