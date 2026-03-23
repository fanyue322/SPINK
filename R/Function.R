#' Remove spatial pattern from gene expression or chromatin accessibility data
#'
#' This function uses thin-plate splines (TPS) to model and remove spatial 
#' dependence from molecular data, creating spatially decorrelated residuals.
#'
#' @param y Numeric vector of expression or accessibility values
#' @param locations Data frame with spatial coordinates (x, y)
#' @param k Dimension of the basis used to represent the smooth term (default: 300)
#' @param method Smoothing parameter estimation method (default: 'GCV.Cp')
#'
#' @return Numeric vector of spatially decorrelated residuals
#' 
#' @import mgcv 
#' @keywords internal
RemoveSpatialPattern<-function(y,location,k=300,method='GCV.Cp')
{
if (length(y) != nrow(location)) {
   stop("Length of y must match number of rows in locations")
}
if (ncol(location) < 2) {
  stop("locations must have at least 2 columns for x and y coordinates")
}
  # Create data frame for GAM
sim_data=data.frame(sx=location[,1],sy=location[,2],X=y)
# Fit thin-plate spline model
# fx=TRUE means fixed degrees of freedom (not penalized)
f_y_hat<-gam(X~s(sx,sy,k=k,fx=TRUE),data=sim_data,method=method)$fitted.values #alternative:method="REML"
 # Calculate residuals (spatially decorrelated)
r_Y<-y-f_y_hat
return(r_Y)
}

#' Preprocess spatial data by removing spatial patterns for all features
#'
#' Applies RemoveSpatialPattern to all genes or peaks in parallel.
#'
#' @param expression.data Matrix of expression or accessibility data (features x spots)
#' @param locations Data frame with spatial coordinates
#' @param num.core Number of cores for parallel processing (default: 1)
#'
#' @return Matrix of spatially decorrelated data with same dimensions as input
#' 
#' @importFrom pbmcapply pbmclapply
#' @importFrom parallel detectCores
#' @keywords internal
SpatialPreprocess<-function(expression.data,location,num.core=1)
{

index=1:nrow(expression.data)

res.exp <- pbmcapply::pbmclapply(index, mc.cores = num.core, function(i){
          gene.expression <- as.numeric(x = expression.data[i, , drop = FALSE])
		  print(paste0("####perform gene expression spatial pattern remove for Gene",i,'####'))
		  tryCatch({suppressWarnings(
		    res<-RemoveSpatialPattern(gene.expression,location))
  }, warning=function(w){ 
		      print(w); return(res);
		    }, error=function(e){ 
		      print(e); return(NULL);
		    }, finally={
		      return(res)
		    }
		  )
		})
res.exp<-do.call(rbind,res.exp)	
rownames(res.exp)=rownames(expression.data)
colnames(res.exp)=colnames(expression.data)
return(res.exp)
}

#' Generate permutation peaks matched on genomic features
#'
#' For each peak, finds 200 background peaks from different chromosomes
#' matched on GC content, total accessibility, and sequence length.
#'
#' @param peak.name Name of the query peak
#' @param all.peaks Vector of all available peak names
#' @param meta.features Data frame with peak metadata (GC.percent, count, sequence.length)
#'
#' @return List of matched background peak names
#' 
#' @importFrom stringr str_split
#' @importFrom Signac MatchRegionStats
#' @keywords internal
GeneratePermutePeak<-function(peak.name,all.peaks,meta.features)
{
features.match <- c("GC.percent", "count", "sequence.length")
gene.chrom=stringr::str_split(peak.name,'-')[[1]][1]
tran.peaks <- all.peaks[!grepl(pattern = paste0("^", gene.chrom), x = all.peaks)]
                trans.peaks <- all.peaks[!grepl(pattern = paste0("^", gene.chrom), x = all.peaks)]
                meta.use <- meta.features[trans.peaks, ]
                pk.use <- meta.features[peak.name, ]
                bg.peaks <- lapply(X = seq_len(length.out = nrow(x = pk.use)), 
                  FUN = function(x) {
                    Signac:::MatchRegionStats(meta.feature = meta.use, 
                      query.feature = pk.use[x, , drop = FALSE], 
                      features.match = features.match, n = 200, 
                      verbose = FALSE)
                  })
return(bg.peaks)
}

#' Fit constrained spline model for enhancer-gene regulation
#'
#' Fits a hierarchical linear model with non-negative constraints on 
#' enhancer effects using cone programming.
#'
#' @param y Gene expression vector (spatially decorrelated)
#' @param x Matrix of chromatin accessibility (peaks x spots)
#' @param z Covariate matrix (optional, default: intercept only)
#' @param nsim Number of simulations for p-value calculation (default: 100)
#'
#' @return List containing:
#'   \item{muhat}{Fitted values}
#'   \item{bh}{Coefficient estimates}
#'   \item{pval}{Overall p-value for regulation}
#'   \item{p_bh}{Peak-specific p-values}
#' 
#' @importFrom coneproj coneB
#' @importFrom stats pchisq pbeta rnorm sd
#' @keywords internal
RunSpatialATAC_gSEM_NonNegative<-function(y,x,z=NULL,nsim=100)
{
require(coneproj)

 xmat=as.matrix(x)
 n = length(y)
 if(is.null(z))
 {
 z=rep(1,n)
 }else{
 z=cbind(rep(1,n),z)
 }
 zmat=as.matrix(z)
 zvec = y
 dsend = xmat
 zsend = zmat
 ans = coneB(zvec, dsend, zsend)
 bigmat=cbind(zsend,dsend)
 edf = ans$df
            face = ans$face
            bh = coef(ans)
			np=1
			    if (any(round(bh[1:np],6) < 0)) {
                pos = (1:np)[which(round(bh[1:np],6) < 0)]
                face = unique(c(pos, face))
            }
			
			 dd = t(bigmat[face, ,drop = FALSE])
			
			muhat = bigmat %*% bh
			
			yhat=muhat                                                    
			pmat=zsend%*%solve(t(zsend)%*%zsend)%*%t(zsend)
			th0=pmat%*%y
			sse0=sum((y-th0)^2)
			sse1=sum((y-yhat)^2)
		bstat=(sse0-sse1)/sse0
		m=ncol(bigmat)
		mdist=1:(m+1)*0
		k0=dim(zsend)[2]
		for(isim in 1:nsim){
			ysim=rnorm(n)
			asim=coneB(ysim,dsend,zsend)
			df0=asim$df-k0
			mdist[df0+1]=mdist[df0+1]+1
		}
		mdist=mdist/nsim
		ps=mdist[1]
		for(d in 1:m){
			ps=ps+pbeta(bstat,d/2,(n-d-k0)/2)*mdist[d+1]
		}
		pval=1-ps
		sig = 1/(n-ans$df)*sum((zvec-bigmat%*%bh)^2)
		logLik<-sum(log(dnorm(x=zvec,mean=bigmat%*%bh,sd=sqrt(sig))))
		if((length(bh)-ncol(zmat))>1)
		{
		p_bh=c()
		for(i in (ncol(zmat)+1):nrow(bh))
		{
		ans = coneB(zvec, as.matrix(dsend[,-(i-ncol(zmat))]), zsend)
		bigmat=cbind(zsend,dsend[,-(i-ncol(zmat))])
		sig = 1/(n-ans$df)*sum((zvec-bigmat%*%ans$coefs)^2)
		logLik2<-sum(log(dnorm(x=zvec,mean=bigmat%*%ans$coefs,sd=sqrt(sig))))
		tmp=2*(logLik-logLik2)
		ps=pchisq(tmp,df=1,lower.tail=F)
		p_bh=c(p_bh,ps)
		}
		}else{
		p_bh=pval
		}
    rslt = list(muhat = muhat,bh = bh,pval=pval,p_bh=p_bh)
    return (rslt)
}

#' Get random permutation peaks for null distribution
#'
#' @param peak.access Vector of peak names to permute
#' @param permute_peak List of pre-computed permutation peaks
#'
#' @return Vector of randomly selected background peaks
#' @keywords internal
GetPermutePeak<-function(peak.access,permute_peak)
{
idx=match(peak.access,names(permute_peak))
permute_peak.use<-c()
for(i in idx)
{
num=length(permute_peak[[i]][[1]])
permute_peak.use<-c(permute_peak.use,permute_peak[[i]][[1]][sample(1:num)[1]])
}

return(permute_peak.use)
}

#' Core function to link peaks to genes
#'
#' Performs gene-level constrained regression and peak-level correlation analysis.
#'
#' @param object Seurat object
#' @param res Preprocessing results list
#' @param n.permute Number of permutations
#' @param cor_method Correlation method (default: qlcMatrix::corSparse)
#'
#' @return Data frame with gene-level and peak-level results
#' 
#' @importFrom pbmcapply pbmclapply
#' @importFrom qlcMatrix corSparse
#' @keywords internal
LinkSpatialPeaks<-function(object,res,n.permute=1000)
{
  expression.data.processed<-res$processed.rna
  peak.data.processed<-res$processed.peak
  permute_peak<-res$permute_peak
  peak_distance_matrix<-res$peak_distance_matrix
  spot.use<-res$spot.use
  genes.use<-res$genes.use
  all.peaks<-res$all.peaks
  num.core<-res$num.core
    cor_method <- qlcMatrix::corSparse
########Gene Level analysis##########
  result.rna=data.frame()
  for(i in 1:length(genes.use))
  {
  print(paste0("processing Gene",genes.use[i]))
  peak.use <- as.logical(x = peak_distance_matrix[, genes.use[i]])
  expression.data.x <- as.numeric(x = expression.data.processed[i,spot.use , drop = FALSE])
  peak.access <- peak.data.processed[ peak.use,spot.use]
  if(class(peak.access)[1]=='numeric')
  {
  peak.access=t(as.matrix(peak.access))
  }else{
  peak.access=(as.matrix(peak.access))
  }
  res<-RunSpatialATAC_gSEM_NonNegative(expression.data.x,t(peak.access),z=NULL)
  res.sl <-pbmcapply::pbmclapply(1:n.permute, mc.cores = num.core, function(x){
  peak.use.null=GetPermutePeak(all.peaks[peak.use],permute_peak)
  peak.access.null<-peak.data.processed[peak.use.null,spot.use]
  if(class(peak.access.null)[1]=='numeric')
  {
  peak.access.null=t(as.matrix(peak.access.null))
  }else{
  peak.access.null=(as.matrix(peak.access.null))
  }
  res<-RunSpatialATAC_gSEM_NonNegative(expression.data.x,t(peak.access.null),z=NULL)
  return(res$pval)
  })
  pnull=do.call(c,res.sl)
  p_permute=length(which(pnull<res$pval))/n.permute
  result.rna=rbind(result.rna,data.frame(gene=genes.use[i],npeaks=nrow(peak.access),pval=res$pval,p.adj=p_permute))
  }
  ########Peak-gene pair Level analysis##########
   res <- pbmcapply::pbmclapply(X = 1:length(genes.use),mc.cores = num.core, FUN = function(i) {
   peak.use <- as.logical(x = peak_distance_matrix[, genes.use[i]])
   expression.data.x <- (x = expression.data.processed[genes.use[i],spot.use , drop = FALSE])
   peak.access <- peak.data.processed[ all.peaks[peak.use],spot.use, drop = FALSE]

   coef.result <- cor_method(X = t(peak.access), Y = t(expression.data.x))
            rownames(x = coef.result) <- rownames(x = peak.access)
            bg.peaks=permute_peak[all.peaks[peak.use]]
			min_lengths <-min(sapply(permute_peak[all.peaks[peak.use]], function(x) length(x[[1]])))
		     bg.access <-peak.data.processed[ unlist(x = bg.peaks),spot.use, drop = FALSE]
                bg.coef <- cor_method(X = t(bg.access), Y = t(expression.data.x))
                rownames(bg.coef) <- rownames(bg.access)
                zscores <- vector(mode = "numeric", length =  nrow(peak.access))
				n_sample=min(200,min_lengths)
                for (j in seq_along(along.with = rownames(peak.access))) {
                  coef.use <- bg.coef[(((j - 1) * n_sample) + 1):(j * n_sample), ]
                  z <- (coef.result[j] - mean(x = coef.use))/sd(x = coef.use)
                  zscores[[j]] <- z
                }
				pval <- pnorm(q = -abs(x = zscores))     
                d=data.frame(gene=rownames(expression.data.processed)[i],peak=	rownames(peak.access),zscore=zscores,pvalue=pval)	
                return(d)
		})
 result.peak=do.call(rbind,res)
 result.peak$pvalue=result.peak$pvalue*2
 result<-merge(result.rna,result.peak,by='gene')
 colnames(result)<-c("gene","npeaks","pval.gene","p.adj.gene","peak","zscore.peak","pval.peak")
 return(result)
}
 
#' Retrieve SPINK linkage results from Seurat object
#'
#' @param object Seurat object after spink_analysis
#' @param domain Character string specifying which domain results to retrieve.
#'   If NULL (default), returns the first available result. If specified,
#'   returns results for that specific domain.
#'
#' @return Data frame with linkage results
#' 
#' @export

GetLinkResult <- function(object, domain = NULL) {
  spink_data <- object@misc$spink
  
  if (is.null(spink_data)) {
    stop("No SPINK data found. Run spink_preprocess first.")
  }
  
  link_names <- grep("^Link_", names(spink_data), value = TRUE)
  
  if (length(link_names) == 0) {
    stop("No SPINK linkage results found. Run spink_analysis first.")
  }
  
  # If domain specified, look for exact match
  if (!is.null(domain)) {
    target_name <- paste0("Link_", domain)
    if (!target_name %in% link_names) {
      stop("Domain '", domain, "' not found. Available domains: ", 
           paste(gsub("^Link_", "", link_names), collapse = ", "))
    }
    result <- spink_data[[target_name]]
  } else {
    # Return first available result
    result <- spink_data[[link_names[1]]]
    message("Returning results for domain: ", gsub("^Link_", "", link_names[1]))
  }
  
  return(result)
}