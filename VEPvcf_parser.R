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
    stringr::str_extract(pattern = '##INFO=\\<(.+)\\>', group = 1)
  if(length(vepInfo) == 0){
    warning('VEP info header not found. Cannot populate vep fields.')
    return(NA)
  }
  vepInfo %<>% stringr::str_split_1(pattern = ',') %>% str_split('=', simplify = T)
  vepDescr <- vepInfo[vepInfo[,1]=="Description",2]
  vepDescr %<>% str_extract(pattern='Format: (.+)\"' , group = 1)
  return(vepDescr %>% stringr::str_split_1(pattern ="\\|") %>% trimws())
}

getSnpEffcolnames_fromMeta <- function(meta){
  snpEff <- meta[grepl(pattern = '##SnpEffVersion=',meta)]
  snpEffInfo <- meta[grepl(pattern = '##INFO=<ID=ANN',meta)] %>% 
    stringr::str_extract(pattern = '##INFO=\\<(.+)\\>', group = 1) 
  if(length(snpEffInfo) == 0){
    warning('snpEff info header not found. Cannot populate snpEff fields.')
    return(NA)
  }
  snpEffInfo %<>% stringr::str_split_1(pattern = ',') %>% str_split('=', simplify = T)
  snpEffDescr <- snpEffInfo[snpEffInfo[,1]=="Description",2]
  snpEffDescr %<>% str_extract(pattern="Functional annotations: \'(.+)\'" , group = 1)
  return(snpEffDescr %>% stringr::str_split_1(pattern ="\\|") %>% trimws())
}

parseAnnoText <- function(x, annColnames, fieldName = 'CSQ'){  # for snpEff, set to 'ANN'
  if(is.na(annColnames[1]) | is.na(x)) {
    return(NA)
  }
  nc=length(annColnames)
  info <- stringr::str_split_1(x, pattern=';') %>% str_split('=', simplify = T)
  filt = which(info[,1]==fieldName)
  if(length(filt) == 0){ # usually just a spanning deletion (ALT= '*')
    return(NA)
  } else if(length(filt)> 1){
    warning('Multiple ',fieldName,' fields found. just using first')
    filt=filt[1]
  }
  ann = info[filt,2] %>% 
    str_split_1(",") %>% 
    str_split_fixed("\\|", n = nc)
  colnames(ann)=annColnames
  return(ann)
}


melt_VEP <- function(veps, 
                     vcf=NULL, 
                     handles=NULL, 
                     outCols=c('Allele','Consequence', 'IMPACT','SYMBOL','Gene'), 
                     onlyGenes = T) { # from a list of VEP_matrices, 'handles' should be the rownames from which they came; these will be added in a new column
  outCols %<>% as.character()
  if(! is.null(vcf)){
    if(length(veps) != nrow(vcf)){
      warning('VCF table does not have same number of rows as length of vep list of matrices. Ignoring')
      vcf=NULL
    }
  }
  
  if(is.null(handles)){ # just use the original's rownames as a handle
    if(is.null(vcf)){
      names(veps) <- as.character(1:length(veps))
    } else {
      names(veps) <- rownames(vcf)
    }
  } else {
    names(veps) <- handles
  }
  NA_filt <- lapply(veps, function(x){!is.null(dim(x))}) %>% unlist()
  if(!is.null(vcf)){
    Spanning_del <- vcf$ALT != '*'
    if(any(Spanning_del != NA_filt)){
      warning('Some entries in VCF have no VEP matrix and are not spanning deletions.')
    }
  }
  veps <- veps[NA_filt]
  
  # use first element in geneColChoices that matches any of the colnaes of the first vep matrix 
  geneColChoices <- c('Gene', 'Gene_Name', 'Gene_ID', 'SYMBOL')  # VEP 
  geneCol <- geneColChoices [ which(geneColChoices %in% colnames(veps[[1]]))[1] ] 
  if(onlyGenes & length(geneCol) > 0){
    GENEFILT = ''  # <-- this will filter out all effects where Gene == ''
    veps %<>% lapply(FUN=function(x){
      if(length(x) == 1 && is.na(x)){
        return(NA)
      } 
      else {
        return(x[x[,geneCol] != GENEFILT, outCols, drop=F])   # <--- Filter out "NO GENE" effects (GENEFILT='' by default)
      }
    }) 
  } else {
    if( length( geneCol ) == 0 & onlyGenes){
      warning('Non-standard field name for gene in annotation fields. Cannot filter out non-gene associated effects.')
    }
    veps %<>% lapply(FUN=function(x){
      if(length(x) == 1 && is.na(x)){
        return(NA)
      } 
      else {
        return(x[,outCols, drop=F])
      }
    })
  }
  long <- do.call(rbind, veps) %>% as.data.frame()
  long$original_rowname <- rep(names(veps), times=unlist(lapply(veps, nrow)))
  return(long)
}


readMultiVcfs <- function(vcfFiles, sampleNames, sampleNameColumn = 'SAMPLEID', filterIn = NULL, parseAnno = 'none'){ #'VEP' or 'SNPEFF'
  nf <- length(vcfFiles)
  if(nf != length(sampleNames)){
    stop('Number of sampleNames not equal to number of files.')
  } 
  allVcfs <- list()
  allSamps <- list()
  allDP <- list()
  allGQ <- list()
  #allAD <- list()
  allGT <- list()
  allAD_ref <- list()
  allAD_alt <- list()
  allInfo <- list()
  allVEPmats <- list()
  allSnpEffmats <- list()
  vepColnames <- list()
  snpEffColnames <- list()
  for(i in 1:nf){
    sn_i <- sampleNames[i]
    cat('\n processing file ', vcfFiles[i], '\n')
    vepList <- vcf2list(fileName = vcfFiles[i], filterIn = filterIn) #, sn = sn_i, sampleNameColumn = sampleNameColumn)
    vepList$variants[[sampleNameColumn]] = sn_i  # <-- add sample name as last column
    allVcfs[[ sn_i ]] <- vepList$variants
    allSamps[[ sn_i ]] <- vepList$samples
    allDP[[ sn_i ]] = vepList$DP
    allGQ[[ sn_i ]] = vepList$GQ
    allGT[[ sn_i ]] = vepList$GT
    #allAD[[ sn_i ]] = vepList$AD
    allAD_ref[[ sn_i ]] = vepList$REF
    allAD_alt[[ sn_i ]] = vepList$ALT
    allInfo[[ sn_i ]]  = vepList$info
    vepColnames[[ sn_i ]] = vepList$vepColnames
    snpEffColnames[[ sn_i ]] = vepList$snpEffColnames
    
    # It's a good idea to check that all vep colnames are the same, since we imported different files
    if(i>1 && length(vepColnames) > 1){
      diffFields <- setdiff(vepList$vepColnames, vepColnames[[ length(vepColnames)-1 ]])
      if(length(diffFields) > 0){
        warning(paste0('VEP field names do not appear to be consistent for entry ',sn_i, ' with previous.  \n Offending names: ', diffFields, collapse = ',' ))
      }
    }
    
    if('VEP' %in% toupper(parseAnno)){
      cat('\n  * Processing VEP fields under CSQ field')
      allVEPmats[[ sn_i ]] <- extractAnn(vcf = vepList$variants,   # TODO extractAnn would be more efficient if could input vepList$INFO$CSQ
                                         annColnames = vepList$vepColnames, fieldName = 'CSQ',
                                         varHandles = vepList$variants$varHandle )
    }
    if ('SNPEFF' %in% toupper( parseAnno )) {
      cat('\n  * Processing SNPEFF fields under ANN field')
      allSnpEffmats[[ sn_i ]] <- extractAnn(vcf = vepList$variants,   # TODO extractAnn would be more efficient if could input vepList$INFO$CSQ
                                         annColnames = vepList$snpEffColnames, fieldName = 'ANN',
                                         varHandles = vepList$variants$varHandle )
      
    }
  }
  cat('\n Finished importing.\n')
  
  # It's a good idea to check that all vep colnames are the same, since we imported different files
  common_cols <- Reduce(intersect, vepColnames)
  for(i in 1:length(vepColnames)){
    diffFields <- vepColnames[[i]] [ ! vepColnames[[i]] %in% common_cols]
    if(length(diffFields) > 0){
      warning('Extra VEP field names in VCF:', vcfFiles[i],
                ' \n Offending names: ', 
                paste0(diffFields, collapse = ','), '\n')
    }
  }
  return(list(allVcfs = allVcfs, 
              allInfo = allInfo,
              allSamps = allSamps, 
              vepColnames = vepColnames, 
              snpEffColnames = snpEffColnames,
              allVEPmats = allVEPmats,
              allSnpEffmats = allSnpEffmats,
              fmats=list(allDP = allDP, allGQ = allGQ, allGT = allGT, allAD_ref = allAD_ref, allAD_alt = allAD_alt)))
}


vcf2list <- function(fileName, filterIn = NULL, formFields = c('DP','GT','GQ','AD')) {  #, sn = NULL, sampleNameColumn = NULL){
  tmp <- vcfR::read.vcfR(fileName, verbose = T)
  output <- list()
  if(! is.null(filterIn)){
    rowFilt <- which(tmp@fix[,'FILTER'] %in% filterIn) 
    if(length(rowFilt) > 0){
      tmp@fix <- tmp@fix[rowFilt,]
      tmp@gt <- tmp@gt[rowFilt,]
    } else {
      warning('No FILTER values passed filter. Ignoring filterIn argument.')
    }
  }
  if( any(grepl(pattern = ',', x = tmp@fix[,'ALT']))){
    cat('\n WARNING: ALT field contains comma-separated alleles on same row.
        \n Multiallelic rows are not tested with VEPvcf_parser.
        \n Use "bcftools norm -m - " (note the minus sign) to split multiallelic sites across rows. \n')
    warning('ALT field contains comma-separated alleles on same row. Wide format VCFs are not tested with VEPvcf_parser.')
  }
  cat ('\n  * Processing info field')
  output[['info']] <- vcfR::extract_info_tidy(tmp)  %>% as.matrix()  %>% apply(., 2, trimws) # trim whitespace
  cat('\n  * Extracting variant fields')
  v_i <- as.data.frame(tmp@fix) #cbind(tmp@fix, tmp@gt)) 
  v_i$varHandle <- apply(v_i[,c('CHROM', 'POS', 'REF', 'ALT')], 1, paste0, collapse='_')
  cat('\n  * Processing variant, samples, DP, GQ, GT fields')
  output[['variants']] <- v_i
  if(nrow(tmp@gt) > 0){
    output[['samples']] <- as.data.frame(tmp@gt)
    colnames(output[['samples']]) <- make.names( colnames(output[['samples']]) )
    # assumes DP, AD, GQ and GT are present, and only 2 comma-separated values in AD (for ref and alt)
    for(ff in formFields){
      output[[ ff ]] <-  vcfR::extract.gt(x = tmp, element = 'AD', IDtoRowNames = F, as.numeric = ff %in% c('DP', 'GQ')) %>% as.matrix()
      if(ff == 'AD'){
        cat('\n  * Processing AD field into separate REF and ALT readcount matrices, assuming only 1 ALT allele per line. Others will be ignored.')
        tmp_AD <- output[['AD']]
        tmp_AD2 <- apply(tmp_AD, 2, function(x){as.numeric(stringr::str_split_fixed(x,',',n = 2))}) %>% as.matrix()
        if(nrow(tmp_AD2) != nrow(tmp@gt) *2){
          cat('\n WARNING: nrow AD =', nrow(tmp_AD2), ' but nrow @gt is ', nrow(tmp@gt), '\n')
        } else {
          output[['REF']] <- tmp_AD2[1:nrow(tmp@gt),,drop=F]
          output[['ALT']]<-  tmp_AD2[-c(1:nrow(tmp@gt)),,drop=F]
        }
      }
    }
  } else {
    cat('\n no @gt object available.')  
  }
  for(nm in names(output)){
    rownames(output[[nm]]) <- v_i$varHandle  # in theory varHandle should be unique per-row since it is the pasting together of CHROM,POS,REF and ALT
  }
  output[['vepColnames']] <- getVEPcolnames_fromMeta(tmp@meta)
  output[['snpEffColnames']] <- getSnpEffcolnames_fromMeta(tmp@meta)
  return(output)
}


extractAnn <- function(vcf, annColnames, varHandles = NULL, fieldName = 'CSQ'){  # fieldname = 'ANN' for snpEff
  allRes <- list()
  for (i in 1:nrow(vcf)){
    if(vcf$ALT[i] == '*' ){  # spanning deletion: normally no CSQ field
      tmp=NA
    } else {
      tmp = parseAnnoText(x = vcf$INFO[i],  fieldName = fieldName,
                         annColnames = annColnames)
      #if(any(is.na(tmp))) { warning( paste0("No CSQ field found for entry number: ", i, "length:", length(tmp)) )}  #<<>>
    }
    allRes[[i]] <- tmp
  }
  if(! is.null(varHandles)){
    names(allRes) <- varHandles # apply(v_i[,c('CHROM', 'POS', 'REF', 'ALT')], 1, paste0, collapse='_')
  }
  return(allRes)  # list of length (nrow(vcf)) containing VEP matrices; each row in each matrix is 1 effect. 
}


