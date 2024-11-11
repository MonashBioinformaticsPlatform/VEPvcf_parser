# This is an example script, meant for importing vcf files that are chunked over different regions of the genome (or separate INDEL and SNP vcfs), 
# but NOT multiple vcfs of same region from different samples, as it assumes none of the variants are shared between vcf files.
# Also it's only tested on one-alt-allele-per-line vcfs, not comma-separated alt alleles

rbind_VEP_VCF_files <- function(vcfFiles, 
                                sampleNames = NULL, ...){      # if not given, will use vcf file name stems as the sampleNames

  # decide what sample short names should be, based on filenames if not specified in sampleNames
  if(!is.null(sampleNames)){
    if(length(sampleNames) != length(vcfFiles)){
      stop('sampleNames, if specified, must have same number of elements as vcfFiles\n')
    }
    filesXsamps = data.frame(file=vcfFiles, samp = sampleNames)
  } else {
    filesXsamps = data.frame(file=vcfFiles)
    filesXsamps$samp = vcfFiles %>%
      base::basename() %>%
      stringr::str_remove('.gz$') %>%
      stringr::str_remove('.vcf$')
  }
  filesXsamps$samp2 <-  make.names(filesXsamps$samp, unique = T)
  
# import vcf files and parse VEP parts of INFO column
  vcfs_imp_l <- readMultiVcfs(vcfFiles = vcfFiles,  
                              sampleNames = filesXsamps$samp2,
                              sampleNameColumn="SAMPLEID",    # SAMPLEID is default colname produced by mutSignatures::importVCFfiles (to maintain compatibility)
                              filterIn = 'PASS',              # Remove rows unless PASS in filter field to reduce size of table. Default = no pre-filter
                              parseAnno = c('VEP','SNPEFF'))  # Assuming that the vcf file contains both VEP and SNPEFF output

  # It's a good idea to check that all vep colnames are the same, since we imported different files
  common_cols <- Reduce(intersect, vcfs_imp_l$vepColnames)
  all_cols <- vcfs_imp_l$vepColnames %>% unlist() %>% unique()
  cat('\n Shared VEP field names: ',paste0(common_cols, sep=','),'\n')
  if(length(all_cols) > length(common_cols)) {
    warning('There were some VEP fieldnames not present in all files.')
  }

  # rbind all the results from different vcf files
  cat('\n rbinding vcf tables / matrices \n')
  vcfs_imp <- do.call(rbind, vcfs_imp_l$allVcfs)
  samps_imp <- do.call(rbind, vcfs_imp_l$allSamps)
  info_imp <-  do.call(rbind, vcfs_imp_l$allInfo)

  # check if order somehow messed up
  if(all(vcfs_imp$SAMPLEID == samps_imp$SAMPLEID) & all(rownames(vcfs_imp) == rownames(samps_imp))){
    vcfs_imp %<>% cbind(., samps_imp) #<--- only works if patient / column names are compatible and proper names
  } else {
    warning('Sample order somehow messed up? Not cbinding vcfs with samples.')
  }

  # rbind all the fmats matrices of sample-specific values
  fmats <- list()
  for(nm in names(vcfs_imp_l$fmats)){
    fmats[[nm]] <- do.call(rbind, vcfs_imp_l$fmats[[nm]])
  }

  # check rownames
  x <- lapply(fmats, rownames) %>% as.data.frame()
  if(!all(apply(x, 1, function(x){length(unique(x)) == 1}))){
    warning('Scrambled rows in fmats.')
  }
  cat('\n Finished rbinding \n')

  # join VEP / snpEff lists from different files
  conCat_veps <- do.call(c, vcfs_imp_l$allVEPmats)
  conCat_snpEffs <- do.call(c, vcfs_imp_l$allSnpEffmats)

  # Melt VEP / snpEff matrics.
  cat('\n Melting list of matrices to long format \n')
  long_vep_effects_mat <- melt_VEP(veps = conCat_veps, vcf = vcfs_imp, onlyGenes = T, outCols = c('Allele','Consequence', 'IMPACT','SYMBOL','Gene'))
  long_snpEff_effects_mat <- melt_VEP(veps = conCat_snpEffs, vcf = vcfs_imp, onlyGenes = T, outCols = c('Allele', 'Annotation', 'Annotation_Impact','Gene_Name'))
  return(list(vcfs_imp = vcfs_imp, fmats = fmats, long_snpEff_effects_mat = long_snpEff_effects_mat, long_vep_effects_mat = long_vep_effects_mat))
}
