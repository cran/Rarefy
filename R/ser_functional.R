
ser_functional<-function(comm, dist_f=NULL, dist_xy, method=c('rao','chao','fun_div'), tau=NA, q=0, comparison=FALSE, fun_div=NULL, args=NULL, verbose=FALSE) {
  
  method <- method[1]
  if (!method %in% c("rao", "chao", "fun_div")) stop("Unavailable method")
  
  if (!inherits(comm, "matrix") && !inherits(comm, "data.frame")) stop("Non convenient comm")
  if (any(comm<0)) stop("Negative value in comm")
  if(suppressWarnings(any(rowSums(comm))==0)) stop("Empty row")
  if(suppressWarnings(any(colSums(comm))==0)) {
    v<-apply(comm, 2, function(col) any(col !=0 ))
    comm<-comm[,v]
  }
    
  if (!inherits(dist_xy, "dist")) stop("Object of class 'dist' expected for dist_xy") 
  dist_xy<-as.matrix(dist_xy) 
  if (any(dist_xy<0)) stop("Negative value in dist_xy") 
  if (nrow(comm) != nrow(dist_xy)) stop("comm and dist_xy don't have the same number of plots")
  if(!is.null(rownames(dist_xy)) && !is.null(rownames(comm)) ) {
    if(any(!rownames(comm)%in%rownames(dist_xy))) stop("comm and dist_xy must have the same names for the plots")
  } else if(!is.null(rownames(dist_xy)) && is.null(rownames(comm))) {
    rownames(comm)<-rownames(dist_xy)
    warning("comm has no row names")
    warning("row names of dist_xy set as row names of comm")
  } 
  
  
  if(is.null(dist_f) && method!='fun_div') stop("dist_f must have a value")
  if(!is.null(dist_f)) {
    if (!inherits(dist_f, "dist")) stop("Object of class 'dist' expected for dist_f")
    if (!is.euclid(sqrt(dist_f))) warning("Squared Euclidean or Euclidean property expected for dist_f")
    dist_f<-as.matrix(dist_f)
    if (ncol(comm) != nrow(dist_f)) stop("comm and dist_f don't have the same number of species")
    if(!is.null(rownames(dist_f)) && !is.null(colnames(comm)) && any(!colnames(comm)%in%rownames(dist_f))) stop("Names in dist_f must be the same as col names in comm") 
    if (any(dist_f<0)) stop("Negative value in dist_f") 
  }
  
  if(method=='rao') {
    r_fin<-array(dim = c(nrow(comm), nrow(comm)))
    com<-sweep(comm,1,rowSums(comm),"/")
    nami<-rownames(comm)
    for(i in 1:nrow(comm)) {
      v<-nami[order(dist_xy[,i])]
      x<-com[v,]
      x<-apply(x,2,cummean)
      r_fin[i,]<-apply(x,1,function(newp) t(as.vector(newp)) %*% dist_f %*% as.vector(newp))
    }
    fin<-colMeans(r_fin,na.rm=TRUE)
    IC_up <- fin + (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(ncol(comm))))
    IC_low <- fin - (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(ncol(comm))))
    df<-data.frame(fin,IC_up,IC_low)
    colnames(df)<-c('Rarefaction','IC_up','IC_low')
    if(comparison) {
      nsr<-rare_Rao(comm,as.dist(dist_f))
      df<-data.frame(df,nsr)
    }
  }
  
  else if(method=='chao') {
    dist_f[which(dist_f>tau,arr.ind = T)] <- tau
    nami <- rownames(comm)
    r_fin<-array(dim = c(nrow(comm), nrow(comm)))
    for(i in 1:nrow(comm)) {
      v<-nami[order(dist_xy[, i])]
      x<-comm[v,]
      x<-apply(x,2,cumsum)
      for(j in 1:nrow(comm)) {
        v1<-as.matrix(x[j,])
        a <- as.vector((1 - dist_f/tau) %*% as.vector(v1) )
        d<-x[j,]
        d <- d[a!=0]
        nplus <- sum(d)
        a<-a[a!=0]
        d <- d/a
        if(q==1){
          r_fin[i,j]<-exp(sum(-d*a/nplus*log(a/nplus)))
        }else{
          r_fin[i,j]<-(sum(d*(a/nplus)^q))^(1 / (1-q))
        }
      } 
    }
    fin<-colMeans(r_fin,na.rm=TRUE)
    IC_up <- fin + (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(nrow(comm))))
    IC_low <- fin - (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(nrow(comm))))
    df<-data.frame(fin,IC_up,IC_low)
    colnames(df)<-c('Rarefaction','IC_up','IC_low')
    if(comparison) {
      r_fin<-array(dim = c(99, nrow(comm)))
      for(i in 1:99) {
        v<-sample(1:nrow(comm),nrow(comm))
        x<-comm[v,]
        x<-apply(x,2,cumsum)
        for(j in 1:nrow(comm)) {
          v1<-as.matrix(x[j,])
          a <- as.vector((1 - dist_f/tau) %*% as.vector(v1) )
          d<-x[j,]
          d <- d[a!=0]
          nplus <- sum(d)
          a<-a[a!=0]
          d <- d/a
          if(q==1){
            r_fin[i,j]<-exp(sum(-d*a/nplus*log(a/nplus)))
          }else{
            r_fin[i,j]<-(sum(d*(a/nplus)^q))^(1 / (1-q))
          }
        } 
      }
      fin<-colMeans(r_fin,na.rm=TRUE)
      IC_up <- fin + (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(nrow(comm))))
      IC_low <- fin - (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(nrow(comm))))
      df1<-data.frame(fin,IC_up,IC_low)
      colnames(df1)<-c('Classic Rarefaction','IC_up','IC_low')
      df<-data.frame(df,df1)
    }
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
        if(exists(v[ch]) && v[ch]!= n[i]) { 
          l[[i]]<-get(v[ch])
          i<-i+1
        }
        else if(v[ch] %in% c('FALSE','TRUE','T','F') && v[ch]!= n[i]) {
          l[[i]]<-as.logical(v[ch])
          i<-i+1
        }
        else if(!is.na(as.numeric(v[ch])) && v[ch]!= n[i]) {
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
            l[[i]]<-gsub("\'|\"",'',v[ch])
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
          l[[i]]<-v[ch]
          i<-i+1
        }
      }
      
      names(l)<-n
      
      nami<-rownames(comm)
      r_fin<-array(dim = c(nrow(comm), nrow(comm)))
      for(i in 1:nrow(comm)) {
        v<-nami[order(dist_xy[, i])]
        x<-comm[v,]
        x<-apply(x,2,cumsum)
        l[[1]]<-x
        r_fin[i,]<-do.call(f,l)
      }
      fin<-colMeans(r_fin,na.rm=TRUE)
      IC_up <- fin + (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(nrow(comm))))
      IC_low <- fin - (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(nrow(comm))))
      df<-data.frame(fin,IC_up,IC_low)
      colnames(df)<-c('Rarefaction','IC_up','IC_low')
      if(comparison) {
        r_fin<-array(dim = c(999,nrow(comm)))
        for(i in 1:999) {
          v<-sample(1:nrow(comm),nrow(comm))
          x<-comm[v,]
          x<-apply(x,2,cumsum)
          l[[1]]<-x
          r_fin[i,]<-do.call(f,l)
        }
        fin<-colMeans(r_fin,na.rm=TRUE)
        IC_up <- fin + (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(nrow(comm))))
        IC_low <- fin - (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(nrow(comm))))
        df1<-data.frame(fin,IC_up,IC_low)
        colnames(df1)<-c('Classic Rarefaction','IC_up','IC_low')
        df<-data.frame(df,df1)
      }
    }
    
    else {
      if(!all(names(args) %in% v)) stop("The arguments must be the ones specified by the function")
      
      ind<-match(NA,unlist(a))
      ind<-match(names(unlist(args)[ind]), names(args))
      
      nami<-rownames(comm)
      r_fin<-array(dim = c(nrow(comm), nrow(comm)))
      for(i in 1:nrow(comm)) {
        v<-nami[order(dist_xy[, i])]
        x<-comm[v,]
        x<-apply(x,2,cumsum)
        args[[ind]]<-x
        r_fin[i,]<-do.call(f,args)
      }
      fin<-colMeans(r_fin,na.rm=TRUE)
      IC_up <- fin + (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(nrow(comm))))
      IC_low <- fin - (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(nrow(comm))))
      df<-data.frame(fin,IC_up,IC_low)
      colnames(df)<-c('Rarefaction','IC_up','IC_low')
    }
    if(comparison) {
      r_fin<-array(dim = c(999, nrow(comm)))
      for(i in 1:999) {
        v<-sample(1:nrow(comm),nrow(comm))
        x<-comm[v,]
        x<-apply(x,2,cumsum)
        args[[ind]]<-x
        r_fin[i,]<-do.call(f,args)
      }
      fin<-colMeans(r_fin,na.rm=TRUE)
      IC_up <- fin + (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(nrow(comm))))
      IC_low <- fin - (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(nrow(comm))))
      df1<-data.frame(fin,IC_up,IC_low)
      colnames(df1)<-c('Classic Rarefaction','IC_up','IC_low')
      df<-data.frame(df,df1)
    }
  }
  
  return(df)
}

