#library(mutSignatures)
library(stringr)
library(magrittr)
library(vcfR)

getVEPcolnames <- function(x,n=-1){
  vraw <- readLines(con = x, n = n)
  return(getVEPcolnames_fromMeta(meta=vraw))
}

getVEPcolnames_fromMeta <- function(meta){
  vep <- meta[grepl(pattern = '##VEP=',meta)]
  vepInfo <- meta[grepl(pattern = '##INFO=<ID=CSQ',meta)] %>% 
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
    #warning('no CSQ field found') # stop <<>>
    return(NA)                 #      <<>>
  } else if(length(filt)> 1){
    warning('Multiple CSQ fields found. just using first')
    filt=filt[1]
  }
  csq = info[filt,2] %>% 
    str_split_1(",") %>% 
    str_split_fixed("\\|", n = nc)
  colnames(csq)=VEPcolnames
  return(csq)
}

importVEPVCFfiles <- function(vcfFiles, 
                            sampleNameColumn="SAMPLEID", # mutSignatures::importVCFfiles automatically creates this column
                            filterIn=NULL, 
                            sampleNames=NULL, ...){ # if not given, will use vcf file name stem
  # usually, set filterIn='PASS'
  
  if(!is.null(sampleNames)){
    filesXsamps = data.frame(file=vcfFiles, samp = sampleNames)
  } else {
    filesXsamps = data.frame(file=vcfFiles)
    filesXsamps$samp = vcfFiles %>%
    base::basename() %>%
    stringr::str_remove('.gz$') %>%
    stringr::str_remove('.vcf$')
  }
  filesXsamps$samp2 <-  make.names(filesXsamps$samp, unique = T)
  
  vcfs_imp_l <- readMultiVcfs(vcfFiles = vcfFiles,  # mutSignatures::importVCFfiles
                          sampleNames = filesXsamps$samp2,
                          sampleNameColumn = sampleNameColumn) #populate vcfs_imp$SAMPLEID with a derivative of the vcf file-name(s)
                                            #sampleNameColumn = sampleNameColumn) # <<>>
  vcfs_imp <- vcfs_imp_l$allVcfs
  #VEPcolnames <- list()
  #for(rw in 1:nrow(filesXsamps)){
  #  VEPcolnames[[ filesXsamps$samp[rw] ]] <- getVEPcolnames( filesXsamps$file[rw] )  # different vcfs might have different VEP colnames
  #}
  VEPcolnames <- vcfs_imp_l$vepColnames
  
  # reduce size for efficiency
  if(! is.null(filterIn)){
    vcfs_imp <- vcfs_imp[vcfs_imp$FILTER %in% filterIn,]  
  }
  # parse VEP
  #info <- vcfR::extract_info_tidy(x)
  allRes=list()
  for (i in 1:nrow(vcfs_imp)){
    samp <- vcfs_imp[i,sampleNameColumn]  # which file did this originate from? Need to know so we can set up corresponding VEP fields
    if (! samp %in% names(VEPcolnames)){
      stop( 'samp not matching any VEP colname entry' )
    }
    if(vcfs_imp$ALT[i] == '*' ){  # spanning deletion: normally no CSQ field
      tmp=NA
    } else {
      tmp = parseVEPInfo(x = vcfs_imp$INFO[i],
                               VEPcolnames=VEPcolnames[[ samp ]])
    
      if(any(is.na(tmp))) { warning( paste0("No CSQ field found for entry number: ", i, "length:", length(tmp)) )}  #<<>>
    }
    allRes[[i]] <- tmp
  }
  vcfs_imp$VEP_matrix <- allRes
  return(list(vcfs_imp = vcfs_imp, fmats = vcfs_imp_l$fmats))
}

melt_VEP <- function(x, handle=NULL, outCols=c('Allele','Consequence', 'IMPACT','SYMBOL','Gene')){
  veps <- x$VEP_matrix 
  if(is.null(handle)){ # just use the original's rownames as a handle
    names(veps) <- rownames(x)
  } else {
    names(veps) <- x[,handle]
  }
  NA_filt <- lapply(x$VEP_matrix, function(x){!is.null(dim(x))}) %>% unlist()
  Spanning_del <- x$ALT != '*'
  if(any(Spanning_del != NA_filt)){
    warning('Some entries in VCF have no VEP matrix and are not spanning deletions.')
  }
  veps <- veps[NA_filt]
  veps %<>% lapply(FUN=function(x){
    if(length(x) == 1 && is.na(x)){
      return(NA)
    } 
    else {
      return(x[x[,'Gene'] != '',outCols, drop=F])
    }
  })  # <--- Filter out "NO GENE" effects
  long <- do.call(rbind, veps) %>% as.data.frame()
  long$original_rowname <- rep(names(veps), times=unlist(lapply(veps, nrow)))
  return(long)
}


readMultiVcfs <- function(vcfFiles, sampleNames, sampleNameColumn = 'SAMPLEID'){
  nf <- length(vcfFiles)
  if(nf != length(sampleNames)){
    stop('Number of sampleNames not equal to number of files.')
  } 
  allVcfs <- list()
  allDP <- list()
  allGQ <- list()
  #allAD <- list()
  allGT <- list()
  allAD_ref <- list()
  allAD_alt <- list()
  vepColnames <- list()
  for(i in 1:nf){
    tmp <- vcfR::read.vcfR(vcfFiles[i], verbose = T)
    allVcfs[[i]] = as.data.frame(cbind(tmp@fix, tmp@gt))
    allVcfs[[i]][[sampleNameColumn]] = sampleNames[i]
    vepColnames[[ sampleNames[i] ]] <- getVEPcolnames_fromMeta(tmp@meta)
    allDP[[ sampleNames[i] ]] <- vcfR::extract.gt(x = tmp, element = 'DP', IDtoRowNames = F, as.numeric = T)
    allGQ[[ sampleNames[i] ]] <- vcfR::extract.gt(x = tmp, element = 'GQ', IDtoRowNames = F, as.numeric = T)
    allGT[[ sampleNames[i] ]] <- vcfR::extract.gt(x = tmp, element = 'GT', IDtoRowNames = F, as.numeric = F)
    tmp_AD <- vcfR::extract.gt(x = tmp, element = 'AD', IDtoRowNames = F, as.numeric = F)
    tmp_AD2 <- apply(tmp_AD, 2, function(x){as.numeric(stringr::str_split_fixed(x,',',n = 2))})
    allAD_ref[[ sampleNames[i] ]] <- tmp_AD2[1:nrow(tmp@gt),]
    allAD_alt[[ sampleNames[i] ]] <-  tmp_AD2[-c(1:nrow(tmp@gt)),]
  }
  allVcfs <- do.call(rbind, allVcfs)
  #allDP <- do.call(rbind, allDP)
  #allGQ <- do.call(rbind, allGQ)
  #allAD <- do.call(rbind, allAD)
  #allGT <- do.call(rbind, allGT)
  return(list(allVcfs=allVcfs, vepColnames=vepColnames, fmats=list(allDP = allDP, allGQ=allGQ, allGT=allGT, allAD_ref = allAD_ref, allAD_alt = allAD_alt)))
}
