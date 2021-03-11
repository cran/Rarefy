## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(fig.width=7,fig.height=5,fig.align = "center",comment = "#",collapse = FALSE)

## ----results='hide', message=FALSE, warning=FALSE-----------------------------
#Required packages
require(Rarefy)
require(ade4)
require(adiv)
require(ape)
require(vegan)
require(phyloregion)
require(raster)

## ----Tab01, echo=FALSE, warning=FALSE, message=FALSE--------------------------
require(kableExtra)
Functions<-c('$\\textbf{rare_alpha}$','','$\\textbf{rare_beta}$','','','','$\\textbf{rare_phylo/ser_phylo}$','','','','$\\textbf{ser_functional}$','')
Metrics<-c('HCDT entropy (HCDT)','Hill numbers ($^{q}D$)','Bray-Curtis dissimilarity ($\\beta_{bray}$)','Cody index ($\\beta_{c}$)','Jaccard similarity coefficient ($\\beta_{j}$)','Whittaker\'s species turnover ($\\beta_{w}$)','Barker\'s weighted phylogenetic diversity (PDw)','Faith’s phylogenetic diversity (PD)','Feature diversity indexes (ƒdiv)','Pavoine’s index (Ia)','Chao’s functional beta-diversity index (FD)','Rao’s quadratic entropy (Q)')
Formula<-c('$$HCDT=\\frac{\\left(1-\\sum_{i=1}^{n}p_i^q\\right)}{\\left(q-1\\right)}$$','$$^{q}\\textrm{D}=\\left( \\sum_{i=1}^{S}p_i^q\\right )^{1/(1-q)}$$','$$\\beta_{bray}=\\frac{\\sum_{i}(x_i-x_j)}{\\sum_{i}x_i+x_j}$$','$$\\beta_c=\\frac{\\left[g\\left(H\\right)+l\\left(H\\right)\\right]}{2}$$','$$\\beta_j=\\frac{a}{\\left(\\alpha_1+\\alpha_2-a\\right)}$$','$$\\beta_w=\\frac{\\gamma}{\\alpha}-1$$','$$PDw=B\\times\\frac{\\sum_{i}^{B}{L_iA_i}}{\\sum_{i}^{B}A_i}$$','$$PD=\\sum_{i\\in B} L_i$$','$$^qƒdiv_{Hill}=\\left [\\sum_{i \\in B} L_{i}(p_{i})^{q}]\\right ]^{\\frac{1}{1-q}}$$ $$^qƒdiv_{HCDT}= \\frac{1-\\sum_{i \\in B} L_{i}(p_{i})^q}{q-1}$$ $$^qƒdiv_{Renyi}= \\frac{1}{1-q}log\\left [ \\sum_{i \\in B} L_{i}(p_{i})^q\\right]$$','$$I_a=\\sum_{K=1}^{N}{\\left(t_K-t_{K-1}\\right)H_{a,K}}$$','$$^{q}\\textrm{FD}(\\Delta(\\tau))=\\left ( \\sum_{i=1}^{S} \\nu_{i}(\\tau)\\left(\\frac{a_i(\\tau)}{n_{+}} \\right )^{(1/1-q)} \\right )$$','$$Q\\left(p_i,D\\right)=\\sum_{k=1}^{S}\\sum_{l=1}^{S}p_{k}p_{l}d_{kl}$$')
Description<-c('HCDT entropy (Harvda & Charvat 1967; Daròczy 1970; Tsallis 1988) is a generalization of the standard coefficient of entropy. $q$ is the parameter that regulates the sensitivity to species abundance and $p_i$ is the relative abundance of species $i$. For $q=0$, the index is equal to species richness minus one, for $q$ tending to 1, it is equivalent to Shannon entropy and for $q=2$ it is equivalent to the Gini-Simpson index.','Hill numbers (Hill 1973) is a class of measures that obeys to the replication principle and integrates species richness and species abundances. The parameter $q$, called \'order\', regulates the sensitivity of the index to the species abundance: with $q=0$ the value of the index corresponds to the species richness, with $q$ tending to 1 the measure tends to the exponential of Shannon index, and with $q=2$ it corresponds to the inverse of Simpson index. $p_i$ is the relative abundance of species $i$.','Bray-Curtis dissimilarity (Bray & Curtis 1957) is a pairwise measure of similarity between plots weighted by the abundances of the species. The accumulation curve is calculated with the mean of pairwise dissimilarities among $N$ plots. $x_i$ and $x_j$ are the abundances of the species $x$ in the plots $i$ and $j$.','Cody index (Cody 1975) is defined as the rate at which species are being replaced in censuses at each point on the habitat gradient and is fixed for samples arranged along gradients of environmental change. $g(H)$ is the number of species gained along the habitat gradient $H$ and $l(H)$ is the number of species lost.','Jaccard dissimilarity coefficient (Jaccard 1901, 1912) is a pairwise measure of dissimilarity between plots. The accumulation curve is calculated with the mean of pairwise dissimilarities among $N$ plots. $a$ is the number of species in common between two plots, and $\\alpha_1$ and $\\alpha_2$ are the values of alpha diversity (species richness) of the 2 plots compared.','Whittaker\'s species turnover (Whittaker 1960, 1972) calculates how many times there is a change in species composition among the plots.$\\gamma$ is the species richness over all plots compared and $\\alpha$ the average species richness within a single plot','Barker\'s weighted phylogenetic diversity ($PD_w$) (Barker 2002) is the abundance weighted Faith\'s $PD$: the number of branches is multiplied by the weighted mean branch length, with weights equal to the average abundance of species sharing that branch. $L_i$ is the branch length of the branch $i$, and $A_i$ is the average abundance of the species sharing the branch $i$. $B$ is the number of branches in the tree.','Faith’s phylogenetic diversity ($PD$) (Faith 1992), is defined as the sum of branch lengths in a phylogenetic tree for the assemblage of species. $L_i$ is the branch length of the branch $i$ and $B$ is the number of branches in the tree.','Feature diversity indexes (ƒdiv) (Pavoine & Ricotta 2019) are Hill numbers and the HCDT and Rényi entropies adapted for the calculation of the phylogenetic diversity replacing the species with the units of the branch length in the phylogenetic tree. $L_i$ is the branch length of the branch $i$, $p_i$ is the relative abundance of the species sharing the branch $i$ and $q$ is the scaling constant that weights the importance of rarity of the species. $B$ is the number of branches in the tree.','Pavoine’s index ($I_a$) (Pavoine et al. 2009) calculates phylogenetic diversity partitioned between evolutionary periods and between plots defined in terms of spatial and time units. Tsallis or HCDT entropy (Harvda & Charvat 1967; Daròczy 1970; Tsallis 1988)  (it measures diversity by regrouping individuals into categories) is computed for each period of the phylogenetic tree, from the number of lineages that descend from the period and from the relative abundances summed within these lineages within the focal community. With  $a=0$, HCDT is the richness (number of species) minus one and $I_a$ is Faith\'s $PD$ minus the height of the phylogenetic tree; with $a$ tending to 1 HCDT is a generalization of the Shannon index while with $a=2$ HCDT is the Simpson index and $I_a$ is Rao\'s QE applied to phylogenetic distances between species. To apply $I_a$, the phylogeny must be ultrametric. $H_{a,K}$ is HCDT entropy of order $a$ applied to the period $K$ and $t_K-t_{K-1}$ is the length of the period $K$.','Chao’s functional beta-diversity index ($FD$) (Chao et al. 2019) quantifies the effective number of equally-distinct functional groups in the considered plots at the distinctiveness $\\tau$ threshold. Any two species with functional distance greater than or equal to  $\\tau$, are treated as functionally equally-distinct and as belonging to different functional groups with distance $\\tau$. For each pair of species with functional distance lower than  $\\tau$ but different from zero, only a proportion of individuals is considered functionally equally-distinct, the other proportion of individuals is considered functionally indistinct. If the pairwise distance is equal to zero, the two species are treated as belonging to the same functional group. After dividing the set of species to form functionally indistinct groups, the contribution of every species is quantified and then the $FD$ of order $q$ is calculated using the Hill number of order $q$. $a_{i}(\\tau)$ is the combined abundance of all functionally-indistinct individuals from species $i$, $v_{i}(\\tau)=n_{i}/a_{i}(\\tau)$ represents the attribute contribution of species $i$ for a threshold level $\\tau$ ($n_{i}$ is the abundance of species $i$), $n_+$ is the total number of individuals in the community and $q$ is the parameter that determines the sensitivity of the measure to the relative abundance of the species.','Rao\'s quadratic entropy (Rao 1982) incorporates both the relative abundances of species and a measure of the pairwise functional distances between species. It expresses the average difference between two randomly selected individuals with replacements. $p=(p1,...,p_k,...,S)$ is the vector of relative abundances of species, $S$ is  the  number  of  species, $\\mathbf{D}=(d_{kl})$ is the  matrix of functional dissimilarities  among  species, and $d_{kl}$ is the  functional dissimilarity between species $k$ and $l$.')
tab<-data.frame(Functions,Metrics,Formula,Description)
knitr::kable(tab) %>%
  kable_styling(c("bordered","condensed"),full_width = F) %>%
  column_spec(4, width_min ="10cm") 

## -----------------------------------------------------------------------------
data("duneFVG") #plot/species matrix
data("duneFVG.xy") #plots geographic coordinates

## -----------------------------------------------------------------------------
dist_sp<-dist(duneFVG.xy$tot.xy)

## -----------------------------------------------------------------------------
ser_rarefaction<-directionalSAC(duneFVG$total,dist_sp)

## ----fig01--------------------------------------------------------------------
plot(1:128,ser_rarefaction$N_Exact,xlab="M",ylab="Species richness",ylim=c(0,71),pch=1)
points(1:128,ser_rarefaction$N_SCR,pch=2)
legend("bottomright",legend=c("Classic Rarefaction","Spatially-explicit Rarefaction"),pch=1:2)

## -----------------------------------------------------------------------------
a<-list(NA,'Shannon')
names(a)<-c('comm','method')


## -----------------------------------------------------------------------------
rare_shannon<-rare_alpha(duneFVG$total,method="fun_div",random=999,fun_div='speciesdiv',args=a,mean=TRUE)
rare_shannon_sp<-rare_alpha(duneFVG$total,dist_sp,method="fun_div",random=999,fun_div='speciesdiv',args=a,mean=TRUE,spatial=TRUE)

## ----fig02--------------------------------------------------------------------
plot(rare_shannon[,1],ylab="Shannon index",xlab="Number of sampling units",type="l", ylim=range(rare_shannon,na.rm=TRUE))
lines(rare_shannon[,2],lty=2)
lines(rare_shannon[,3],lty=2)
lines(rare_shannon_sp[,1],col=4)
lines(rare_shannon_sp[,2],lty=2,col=4)
lines(rare_shannon_sp[,3],lty=2,col=4)
legend("bottomright",legend=c("Non spatially-explicit Rarefaction","Spatially-explicit Rarefaction"),lty=1,col=c(1,4))

## -----------------------------------------------------------------------------
data(mite)
data(mite.env)
comm_matrix<-mite 

## -----------------------------------------------------------------------------
beta_directional<-directionalSAC(comm_matrix,mite.env$SubsDens)

## ----fig03--------------------------------------------------------------------
plot(1:70,beta_directional$Beta_M,xlab="M",ylab="Beta diversity",ylim=range(c(beta_directional$Beta_M_dir,beta_directional$Beta_M)))
points(1:70,beta_directional$Beta_M_dir,pch=2)
legend("bottomright",legend=c("Non-directional beta","Directional beta"),pch=1:2)

## -----------------------------------------------------------------------------
data(africa)

## -----------------------------------------------------------------------------
#data are in wgs84, reference system was changed to a projected one (LAEA)
Poli_LAEA<-spTransform(africa$polys,CRSobj='+proj=laea +lat_0=52 +lon_0=10 +x_0=4321000 +y_0=3210000 +ellps=GRS80 +units=m +no_defs')
#we need to fix rownames in order to make them comparable to Poli_LAEA
rownames(africa$comm)<-sapply(rownames(africa$comm),function(x){gsub(pattern = "v",replacement = "", x)})
dista<-dist(coordinates(Poli_LAEA[1:50,]))

## -----------------------------------------------------------------------------
raref<-ser_phylo(as.matrix(africa$comm[1:50,colSums(africa$comm)>0]),africa$phylo,dista,comparison=TRUE) 


## ----fig06--------------------------------------------------------------------
plot(raref[,1],ylab="Faith Index",xlab="Number of sampling units",type="l",ylim=range(raref,na.rm=TRUE))
lines(raref[,2],lty=2)
lines(raref[,3],lty=2)
lines(raref[,4],col=2)
lines(raref[,5],lty=2,col=2)
lines(raref[,6],lty=2,col=2)
legend("bottomright",legend=c("spatially-explicit phylogenetic rarefaction","classic phylogenetic rarefaction"),lty=1,col=1:2)

## -----------------------------------------------------------------------------
data(duneFVG.tr8) #species functional traits
tr8_N<-duneFVG.tr8$traits.nat[,c(1,3,4)]
tr8_D<-data.frame(duneFVG.tr8$traits.nat[,2])
tr8_Q<-duneFVG.tr8$traits.nat[,5:15]
tr8dist<-dist.ktab(ktab.list.df(list(tr8_N,tr8_D,tr8_Q)),type=c('N','D','Q'))

## ----warning=FALSE, message=FALSE---------------------------------------------
rareperm<-rao_permuted(duneFVG$alien,tr8dist)

## ----warning=FALSE, message=FALSE---------------------------------------------
tr8n_N<-duneFVG.tr8$traits.nat[,c(1,3,4)]
tr8n_D<-data.frame(duneFVG.tr8$traits.nat[,2])
tr8n_Q<-duneFVG.tr8$traits.nat[,5:15]
tr8a_N<-duneFVG.tr8$traits.ali[,c(1,3,4)]
tr8a_D<-data.frame(duneFVG.tr8$traits.ali[,2])
tr8a_Q<-duneFVG.tr8$traits.ali[,5:15]
tr8ndist<-dist.ktab(ktab.list.df(list(tr8n_N,tr8n_D,tr8n_Q)),type=c('N','D','Q'))
tr8adist<-dist.ktab(ktab.list.df(list(tr8a_N,tr8a_D,tr8a_Q)),type=c('N','D','Q'))
raren<-rare_Rao(duneFVG$native,tr8ndist)
rarea<-rare_Rao(duneFVG$alien,tr8adist)

## ----fig08--------------------------------------------------------------------
plot(raren[,1], ylab="Rao QE",xlab="Number of sampling units",type="l",ylim=range(raren))
lines(raren[,2],lty=2)
lines(raren[,3],lty=2)
lines(rareperm[,1],col=2)
lines(rareperm[,2],lty=2,col=2)
lines(rareperm[,3],lty=2,col=2)
lines(rarea[,1],col=4)
lines(rarea[,2],lty=2,col=4)
lines(rarea[,3],lty=2,col=4)
legend("bottomright", legend=c("Native species Functional Rarefaction","Standardized Functional Rarefaction","Alien species Functional Rarefaction"),lty=1,col=c(1,2,4))

