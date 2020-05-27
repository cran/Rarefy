
rare_phylo<-function(comm,tree=NULL,method=c("faith", "barker", "Ia", "hill", "tsallis", "renyi","fun_div"),exp=0,resampling=99,fun_div=NULL,args=NULL,verbose=FALSE) {
  
  method <- method[1]
  if(!method%in%c("faith","barker","Ia","tsallis", "hill", "renyi","fun_div")) stop("Unavailable method")
  
  if(is.null(colnames(comm))) stop("comm must have names for columns") 
  if (!inherits(comm, "matrix") && !inherits(comm, "data.frame")) stop("Non convenient comm")
  if (any(comm<0)) stop("Negative value in comm")
  if(suppressWarnings(any(rowSums(comm))==0)) stop("Empty row")
  if(suppressWarnings(any(colSums(comm))==0)) {
    v<-apply(comm, 2, function(col) any(col !=0 ))
    comm<-comm[,v]
  }
  
  if(is.null(tree) && method!='fun_div') stop("tree must have a value")
  if(!is.null(tree)) {
    if(!inherits(tree, "phylo") && !inherits(tree, "phylo4")) stop("tree must be of class phylo or phylo4") 
    if(inherits(tree, "phylo4")) suppressWarnings(tree<-as(tree, "phylo"))
    if(any(!colnames(comm)%in%tree$tip.label)) stop("comm contains tip names that are not available in tree")
    if(is.null(tree$edge.length)) stop("edge lengths are required for tree") 
  }
  
  if(method=='faith') {
    r_fin<-array(dim = c(resampling, nrow(comm)))
    for(i in 1:resampling) {
      v <- sample(1:nrow(comm),nrow(comm))
      x<-comm[v,]
      x<-apply(x,2,cumsum)
      sub_tree<-apply(x,1,function(x) treedata(tree,x[x>0],warnings = FALSE)$phy)
      r_fin[i,]<-unlist(lapply(sub_tree,function(x) sum(x$edge.length)))
    }
    rare<-colMeans(r_fin)
  }
  
  else if(method=='barker') {
    comm<-sweep(comm,2,colSums(comm),"/")
    r_fin<-array(dim = c(resampling, nrow(comm)))
    
      get.leaves<-function(x,st){
        leaves.node<-tips(st,x[2])
      }
    
    for(i in 1:resampling) {
      v <- sample(1:nrow(comm),nrow(comm))
      x<-comm[v,]
      x<-apply(x,2,cumsum)
      sub_tree<-apply(x,1,function(x) treedata(tree,x[x>0],warnings = FALSE)$phy)
      lv<-lapply(sub_tree,function(y) apply(y$edge,1,function(z) get.leaves(z,y)))
      for(j in 1:nrow(x)) {
        v1<-x[j,]
        rel<-unlist(lapply(lv[[j]], function(l) mean(v1[l])))
        r_fin[i,j]<-nrow(sub_tree[[j]]$edge) * ((sum(sub_tree[[j]]$edge.length * rel)) / sum(rel))
      }
    }
    rare<-colMeans(r_fin)
  }
  
  else if(method=='Ia') {
    r_fin<-array(dim = c(resampling, nrow(comm)))
    for(i in 1:resampling) {
      v <- sample(1:nrow(comm),nrow(comm))
      x<-comm[v,]
      x<-apply(x,2,cumsum)
      r_fin[i,]<-pIa(tree,x,exponent = exp)[[1]]
    }
    rare<-colMeans(r_fin)
  }
  
  else if(method=='hill' || method=='tsallis' || method=='renyi') {
    r_fin<-array(dim = c(resampling, nrow(comm)))
    for(i in 1:resampling) {
      v <- sample(1:nrow(comm),nrow(comm))
      x<-comm[v,]
      x<-apply(x,2,cumsum)
      r_fin[i,]<-evodivparam(tree,x,method = method, q=exp)
    }
    rare<-colMeans(r_fin)
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
          l[[i]]<-v[ch]
          i<-i+1
        }
      }
      
      names(l)<-n
      
      r_fin<-array(dim = c(resampling, nrow(comm)))
      for(i in 1:resampling) {
        v<-sample(1:nrow(comm),nrow(comm))
        x<-comm[v,]
        x<-apply(x,2,cumsum)
        l[[1]]<-x
        r_fin[i,]<-do.call(f,l)
      }
      rare<-colMeans(r_fin)
    }
    
    else {
      if(!all(names(args) %in% v)) stop("The arguments must be the ones specified by the function you chose")
      
      ind<-match(NA,unlist(args))
      ind<-match(names(unlist(args)[ind]), names(args))
      
      r_fin<-array(dim = c(resampling, nrow(comm)))
      for(i in 1:resampling) {
        v<-sample(1:nrow(comm),nrow(comm))
        x<-comm[v,]
        x<-as.data.frame(lapply(x,cumsum))
        args[[ind]]<-x
        r_fin[i,]<-do.call(f,args)
      }
      rare<-colMeans(r_fin)
    }
  }
  
  IC_plus <- rare + (1.96*(sd(r_fin)/sqrt(resampling)))
  IC_neg <- rare - (1.96*(sd(r_fin)/sqrt(resampling)))
  df<-data.frame(rare,IC_plus,IC_neg)
  colnames(df)<-c('Rarefaction','IC_plus','IC_neg')
  
  return(df)
  
}

