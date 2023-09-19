library(mutSignatures)
library(stringr)
library(magrittr)

getVEPcolnames <- function(x,n=-1){
  vraw <- readLines(con = x, n = n)
  vep <- vraw[grepl(pattern = '##VEP=',vraw)]
  vepInfo <- vraw[grepl(pattern = '##INFO=<ID=CSQ',vraw)] %>% 
    stringr::str_extract(pattern = '##INFO=\\<(.+)\\>', group = 1) %>%
    stringr::str_split_1(pattern = ',') %>% str_split('=', simplify = T)
  vepDescr <- vepInfo[vepInfo[,1]=="Description",2]
  vepDescr %<>% str_extract(pattern='Format: (.+)\"' , group = 1)
  return(vepDescr %>% stringr::str_split_1(pattern ="\\|"))
}

parseVEPInfo <- function(x, VEPcolnames){
  if(is.na(x)) {
    return(NA)
  }
  nc=length(VEPcolnames)
  info <- stringr::str_split_1(x, pattern=';') %>% str_split('=', simplify = T)
  filt = which(info[,1]=='CSQ')
  if(length(filt) == 0){
    stop('no CSQ field found')
  } else if(length(filt)> 1){
    warn('Multiple CSQ fields found. just using first')
    filt=filt[1]
  }
  csq = info[filt,2] %>% 
    str_split_1(",") %>% 
    str_split_fixed("\\|", n = nc)
  colnames(csq)=VEPcolnames
  return(csq)
}

importVEPVCFfiles <- function(vcfFiles, 
                            sampleNameColumn="SAMPLE_ID", 
                            filterIn=NULL, 
                            sampleNames=NULL, ...){
  # usually, set filterIn='PASS'
  filesXsamps = data.frame(file=vcfFiles)
  if(!is.null(sampleNames)){
    filesXsamps$samp = sampleNames
  } else {
  filesXsamps$samp = vcfFiles %>%
    base::basename() %>%
    stringr::str_remove('.gz$') %>%
    stringr::str_remove('.vcf$')
  }
  VEPcolnames <- list()
  for(rw in 1:nrow(filesXsamps)){
    VEPcolnames[[ filesXsamps$samp[rw] ]] <- getVEPcolnames( filesXsamps$file[rw] )
  }
  vcfs_imp <- mutSignatures::importVCFfiles(vcfFiles = vcfFiles,
                                            sampleNames = filesXsamps$samp,
                                            sampleNameColumn = sampleNameColumn)
  # reduce size for efficiency
  if(! is.null(filterIn)){
    vcfs_imp<-vcfs_imp[vcfs_imp$FILTER %in% filterIn,]  
  }
  # parse VEP
  allRes=list()
  for (i in 1:nrow(vcfs_imp)){
    samp <- vcfs_imp[i,sampleNameColumn]
    allRes[[i]] = parseVEPInfo(x = vcfs_imp$INFO[i],
                               VEPcolnames=VEPcolnames[[ samp ]])
  }
  vcfs_imp$VEP_matrix <- allRes
  return(vcfs_imp)
}

melt_VEP <- function(x, handle=NULL, outCols=c('Allele','Consequence', 'IMPACT','SYMBOL','Gene')){
  veps <- x$VEP_matrix 
  if(is.null(handle)){ # just use the original's rownames as a handle
    names(veps) <- rownames(x)
  } else {
    names(veps) <- x[,handle]
  }
  veps %<>% lapply(FUN=function(x){
    x[x[,'Gene'] != '',outCols, drop=F]
  })  # <--- Filter out "NO GENE" effects
  long <- do.call(rbind, veps) %>% as.data.frame()
  long$original_handle <- rep(names(veps), times=unlist(lapply(veps, nrow)))
  return(long)
}
