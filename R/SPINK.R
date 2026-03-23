#####################################################################
# Package: TDEseq
# Version: 0.0.1
# Modified: 2025-07-28 15:16:02
# Title :  Detecting temporal gene expression changes in the developmental stages of single-cell RNA sequencing studies 
# Authors: Yue Fan and Shiquan Sun
# Contacts: xafanyue@xjtu.edu.cn
#          Xi'an Jiatong University, Department of Biostatistics
######################################################################

#'
#' Preprocess spatial multi-ome data for SPINK analysis
#'
#' This function prepares Seurat objects containing spatial RNA and ATAC data
#' for enhancer-gene linkage analysis. It performs spatial decorrelation, and generates permutation peaks.
#'
#' @param object A Seurat object with RNA and ATAC assays
#' @param group.by Column name in meta.data containing domain annotations
#' @param domain Specific domain to analyze
#' @param refGenome Reference genome: 'hg38', 'hg19', or 'mm10' (default: 'hg38')
#' @param distance Maximum distance for peak-gene pairs (default: 5e5)
#' @param rna.assay Name of RNA assay (default: 'Spatial')
#' @param peak.assay Name of ATAC assay (default: 'peaks')
#' @param min.pct.rna Minimum fraction of spots expressing gene (default: 0.1)
#' @param min.pct.peak Minimum fraction of spots with peak (default: 0.05)
#' @param num.core Number of cores for parallel processing (default: 10)
#'
#' @return Seurat object with preprocessing results stored in misc$spink
#' 
#' @import Seurat
#' @import Signac
#' @importFrom GenomeInfoDb seqinfo seqnames seqlevels seqlevels<- genome genome<-
#' @importFrom pbmcapply pbmclapply
#' @export

spink_preprocess <- function(object,
                  group.by=NULL,
				  domain=NULL,
                  refGenome='hg38',
                  distance = 5e+05,
                  rna.assay='Spatial',
                  peak.assay='peaks',
				  min.pct.rna=0.1,
				  min.pct.peak=0.05,
	              num.core=10				  ) {
				  
#############################Data processing##########################				  

#####Annotate Seurat Object#####				  
if(refGenome=='hg38')
{
if (!require('EnsDb.Hsapiens.v86')) 
{
BiocManager::install('EnsDb.Hsapiens.v86')
}
require(EnsDb.Hsapiens.v86)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "hg38"
}else if(refGenome=='hg19'){
if (!require('EnsDb.Hsapiens.v75')) 
{
BiocManager::install('EnsDb.Hsapiens.v75')
}
require(EnsDb.Hsapiens.v75)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "h19"
}else if(refGenome=='mm10'){
if (!require('EnsDb.Mmusculus.v79')) 
{
BiocManager::install('EnsDb.Mmusculus.v79')
}
require(EnsDb.Mmusculus.v79)
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevels(annotations) <- paste0('chr', seqlevels(annotations))
genome(annotations) <- "mm10"
}	

DefaultAssay(object)<-peak.assay
Annotation(object) <- annotations


if(refGenome=='hg38')
{
if (!require('BSgenome.Hsapiens.UCSC.hg38')) 
{
BiocManager::install('BSgenome.Hsapiens.UCSC.hg38')
}
require(BSgenome.Hsapiens.UCSC.hg38)
genome=BSgenome.Hsapiens.UCSC.hg38
}else if(refGenome=='hg19'){
if (!require('BSgenome.Hsapiens.UCSC.hg19')) 
{
BiocManager::install('BSgenome.Hsapiens.UCSC.hg19')
}
require(BSgenome.Hsapiens.UCSC.hg19)
genome=BSgenome.Hsapiens.UCSC.hg19
}else if(refGenome=='mm10'){
if (!require('BSgenome.Mmusculus.UCSC.mm10')) 
{
BiocManager::install('BSgenome.Mmusculus.UCSC.mm10')
}
require(BSgenome.Mmusculus.UCSC.mm10)
genome=BSgenome.Mmusculus.UCSC.mm10
}	

annot<-Annotation(object[[peak.assay]])
gene.coords <- Signac:::CollapseToLongestTranscript(ranges = annot)

Seurat::DefaultAssay(object)<-peak.assay
object<-RunTFIDF(object)
object<-NormalizeData(object,assay=rna.assay)
peak.data <- Seurat::GetAssayData(object, assay = peak.assay, slot = 'data')
expression.data <- Seurat::GetAssayData(object, assay = rna.assay, slot = 'data')

spot.use=which(object@meta.data[,group.by]==domain)

peakcounts <- rowSums(x = peak.data[,spot.use] > 0)
genecounts <- rowSums(x = expression.data[,spot.use] > 0)
min.cells.peak=round(ncol(peak.data[,spot.use]) * min.pct.peak)
min.cells.rna=round(ncol(expression.data[,spot.use]) * min.pct.rna)
peaks.keep <- peakcounts > (min.cells.peak)
genes.keep <- genecounts > min.cells.rna

peak.data <- peak.data[peaks.keep, ]
expression.data <- expression.data[genes.keep, , drop = FALSE]
    
genes <- rownames(x = expression.data)
gene.coords.use <- gene.coords[gene.coords$gene_name %in% genes,]
peaks <- granges(object[[peak.assay]])
peaks <- peaks[peaks.keep]
standard_chr <- seqnames(genome)
peaks <- peaks[seqnames(peaks) %in% standard_chr]
seqlevels(peaks) <- standard_chr[standard_chr %in% seqlevels(peaks)]


  peak_distance_matrix <- Signac:::DistanceToTSS(
    peaks = peaks,
    genes = gene.coords.use,
    distance = distance
  )
genes.use <- colnames(x = peak_distance_matrix)
all.peaks <- rownames(x = peak_distance_matrix)	

genes.use=genes.use[colSums(peak_distance_matrix)!=0]
peak_distance_matrix=peak_distance_matrix[,colSums(peak_distance_matrix)!=0]
all.peaks=all.peaks[rowSums(peak_distance_matrix)!=0]
peak_distance_matrix=peak_distance_matrix[rowSums(peak_distance_matrix)!=0,]
gene.coords.use <- gene.coords.use[gene.coords.use$gene_name %in% genes.use,]

locations=Seurat::GetTissueCoordinates(object)
expression.data=expression.data[genes.use,]
peak.data=peak.data[all.peaks,]

expression.data.processed<-SpatialPreprocess(expression.data,locations,num.core=num.core)
peak.data.processed<-SpatialPreprocess(peak.data,locations,num.core=num.core)

meta.features <- GetAssayData(object = object, assay = peak.assay, slot = "meta.features")
data.use <- GetAssayData(object = object[[peak.assay]], slot = "counts")
hvf.info <- FindTopFeatures(object = data.use, verbose = FALSE)
hvf.info <- hvf.info[rownames(meta.features), , drop = FALSE]
meta.features <- cbind(meta.features, hvf.info)


all.peaks <- rownames(peak.data.processed)
permute_peak<-list()
for(i in 1:length(all.peaks))
{
peak.name=all.peaks[i]
res<-GeneratePermutePeak(peak.name,all.peaks,meta.features)
permute_peak[[i]]<-res
}
names(permute_peak)=all.peaks

  permute_peak <-pbmcapply::pbmclapply(1:length(all.peaks), mc.cores = num.core, function(x){
  peak.name=all.peaks[x]
  res<-GeneratePermutePeak(peak.name,all.peaks,meta.features)
  return(res)
  })
  names(permute_peak)=all.peaks

res<-list(processed.rna=expression.data.processed,processed.peak=peak.data.processed,
permute_peak=permute_peak,peak_distance_matrix=peak_distance_matrix,
genes.use=genes.use,all.peaks=all.peaks,
distance=distance,spot.use=spot.use,
num.core=num.core)

object@misc$spink<-res

return(object)
}

#'
#' Perform SPINK analysis to identify enhancer-gene regulatory links
#'
#' This function performs gene-level and peak-level statistical testing
#' to identify significant enhancer-gene regulatory interactions.
#'
#' @param object A preprocessed Seurat object from spink_preprocess
#' @param link.gene.thr P-value threshold for gene-level significance (default: 0.1)
#' @param n.permute Number of permutations for empirical p-value (default: 1000)
#'
#' @return Seurat object with results stored in misc$spink$Link_domain
#' 
#' @import Seurat
#' @importFrom pbmcapply pbmclapply
#' @export

spink_analysis <- function(object,link.gene.thr=0.1,n.permute=1000) {

res<-slot(object,"misc")
res<-res[['spink']]

result<-LinkSpatialPeaks(object,res)
result$pval.peak[which(result$zscore.peak<=0)]=1
result$pval.peak[which(result$pval.gene>gene.thr)]=1
result=result[order(result$pval.peak),]
object@misc$spink[[paste0("Link_",domain)]]=result
return(object)
}

#########################################
#             CODE END                  #
#########################################

