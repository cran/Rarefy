
rare_alpha<-function(comm,dist_xy=NULL,method=c("HCDT","hill","fun_div"),q=0,random=99,fun_div=NULL,args=NULL,verbose=FALSE,spatial=FALSE,mean=FALSE) {
  
  method <- method[1]
  if(!method%in%c("HCDT","hill","fun_div")) stop("Unavailable method")
  
  if (!inherits(comm, "matrix") && !inherits(comm, "data.frame")) stop("Non convenient comm")
  if (any(comm<0)) stop("Negative value in comm")
  if(suppressWarnings(any(rowSums(comm))==0)) stop("Empty row")
  if(suppressWarnings(any(colSums(comm))==0)) {
    v<-apply(comm, 2, function(col) any(col !=0 ))
    comm<-comm[,v]
  }
  
  if(spatial && is.null(dist_xy)) stop("Object of class 'dist' expected for dist_xy")
  if(!is.null(dist_xy)){
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
  }
  
  if(method=='HCDT') {
    r_fin<-array(dim = c(ifelse(spatial,nrow(comm),random),nrow(comm)))
    nami<-rownames(comm)
    for(i in 1:ifelse(spatial,nrow(comm),random)) {
      if(spatial) v<-nami[order(dist_xy[,i])]
      else v<-sample(1:nrow(comm),nrow(comm))
      x<-comm[v,]
      x<-apply(x,2,cumsum)
      x<-sweep(x,1,rowSums(x),"/")
      if(q==1) r_fin[i,]<-apply(x,1, function(x) -sum(x[x>0]*log(x[x>0])))
      else r_fin[i,]<-apply(x,1, function(x) (1-sum(x[x>0]**q))/(q-1))
    }
  }
  
  else if(method=='hill') {
    r_fin<-array(dim = c(ifelse(spatial,nrow(comm),random),nrow(comm)))
    nami<-rownames(comm)
    for(i in 1:ifelse(spatial,nrow(comm),random)) {
      if(spatial) v <- nami[order(dist_xy[,i])]
      else v <- sample(1:nrow(comm),nrow(comm))
      x<-comm[v,]
      x<-sweep(x,1,rowSums(x),"/")
      x<-apply(x,2,cummean) 
      if(q==1) r_fin[i,]<-apply(x,1, function(x) exp(-sum(x[x>0]*log(x[x>0]))))
      else r_fin[i,]<-apply(x,1, function(x) sum((x[x>0]**q))**(1/(1-q)))
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
      
      r_fin<-array(dim = c(ifelse(spatial,nrow(comm),random),nrow(comm)))
      nami<-rownames(comm)
      for(i in 1:ifelse(spatial,nrow(comm),random)) {
        if(spatial) v<-nami[order(dist_xy[,i])]
        else v<-sample(1:nrow(comm),nrow(comm))
        x<-comm[v,]
        if(mean) { 
          x<-sweep(x,1,rowSums(x),'/')
          x<-apply(x,2,cummean) 
        }
        else x<-apply(x,2,cumsum)
        l[[1]]<-x
        r_fin[i,]<-do.call(f,l)
      }
    }
    
    else {
      if(!all(names(args) %in% v)) stop("The arguments must be the ones specified by the function")
      
      ind<-match(NA,unlist(args))
      ind<-match(names(unlist(args)[ind]), names(args))
      
      r_fin<-array(dim = c(ifelse(spatial,nrow(comm),random),nrow(comm)))
      nami<-rownames(comm)
      for(i in 1:ifelse(spatial,nrow(comm),random)) {
        if(spatial) v<-nami[order(dist_xy[,i])]
        else v<-sample(1:nrow(comm),nrow(comm))
        x<-comm[v,]
        if(mean) { 
          x<-sweep(x,1,rowSums(x),'/')
          x<-apply(x,2,cummean) 
        }
        else x<-apply(x,2,cumsum)
        args[[ind]]<-x
        r_fin[i,]<-do.call(f,args)
      }
    }
  }
  
  rare<-colMeans(r_fin,na.rm = TRUE)
  IC_up <- rare + (1.96*(sd(r_fin,na.rm = TRUE)/sqrt(ifelse(spatial,nrow(comm),random))))
  IC_low <- rare - (1.96*(sd(r_fin,na.rm=TRUE)/sqrt(ifelse(spatial,nrow(comm),random))))
  df<-data.frame(rare,IC_up,IC_low)
  colnames(df)<-c('Rarefaction','IC_up','IC_low')
  
  return(df)
  
}
