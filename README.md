**Description**

Some simple functions to facilitate extracting SnpEff and VEP fields from VCF files.

Many of these are wrappers for vcfR bioconductor package functions with some additional string parsing. Please credit vcfR.

----

**packages required:**

vcfR

stringr

magrittr

----

**Useage example:**

### Get a vcf file parsed into a list of several results matrices

```
library(stringr)
library(magrittr)
library(vcfR)
source('VEPvcf_parser.R')

# process a vcf file into several different tables:
vcf_in <- vcf2list(fileName = "variants.vcf", filterIn = 'PASS')
```

Most elements in the output list ('vcf_in' above) are tables/matrices having one row per variant.
The names of the elements are: `variants`, `samples`, and (if default `formFields` argument is used): `DP`, `GQ`, `GT`, `REF` and `ALT`.

Also, there should be VEP field names (if present in the file) in `vcf_in$vepColnames` and / or snpEff field names (if present in vcf) in `vcf_in$snpEffColnames` 

### Parse VEP or SnpEff fields:

From the above list, can now parse the VEP fields into a list of matrices. Each matrix corresponds to one variant; each row in the matrix is one **effect** of the variant.

You will also need the parsed vcf from `vcf2list` (here, 'vcf_in').

```
vepMatrices <- extractAnn(vcf = vcf_in$variants, annColnames  = vcf_in$vepColnames, varHandles = vcf_in$variants$varHandle , fieldName = 'CSQ')
snpEffMatrices <- extractAnn(vcf = vcf_in$variants, annColnames  = vcf_in$snpEffColnames, varHandles = vcf_in$variants$varHandle , fieldName = 'ANN')
```

To convert the list of matrices to one long format table (ignoring effects that are not explicitly linked to a gene) use `melt_VEP`:
```
long_effects_mat <- melt_VEP(veps = vepMatrices, vcf = vcf_in, onlyGenes = T, outCols = c('Allele','Consequence', 'IMPACT','SYMBOL','Gene'))
```
This should work for snpEff matrices too **provided** that the nominated `outCols` are fieldnames that are present in the snpEff matrices.

Available fields should be present in `vcf_in$vepColnames` and `vcf_in$snpEffColnames`.

<br>

The output after melting is a matrix, one line per **variant / affected gene** combination, with the columns from the original that are nominated in the `outCols` argument. Therefore, the number of rows will be greater than in the original VCF file.

### Limitations:
This is only tested on one-alt-allele-per-line vcfs, not comma-separated alt alleles.
If you have multiallelic entries on single rows, use "bcftools norm -m - " (note the minus sign) to split multiallelic sites over several rows.
Also it expects some format fields that might not be present in some vcfs ('DP', 'GQ', 'GT' and 'AD'), if these are absent then you will need to specify `formFields` argument with vaules that can be accepted by `vcfR::extract.gt()`. 

### Do lots of files at once:

If vcf files are from non-overlapping genomic regions, with the same set of samples (e.g. results from chunking a file over genomic regions), can use `readMultiVcfs`.

See the below helper function for useage.


```
# This is an example script, meant for importing vcf files that are chunked over different regions of the genome (or separate INDEL and SNP vcfs), 
# but NOT multiple vcfs of same region from different samples, as it assumes none of the variants are shared between vcf files.
# Also it's only tested on one-alt-allele-per-line vcfs, not comma-separated alt alleles

rbind_VEP_VCF_files <- function(vcfFiles, 
                                sampleNameColumn="SAMPLEID", # SAMPLEID is default colname produced by mutSignatures::importVCFfiles (to maintain compatibility)
                                filterIn=NULL,               # usually, set to 'PASS'
                                sampleNames=NULL, ...){      # if not given, will use vcf file name stems as the sampleNames
  
  # decide what sample names should be, based on filenames if not specified in sampleNames
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
                              sampleNameColumn = sampleNameColumn,
                              filterIn = 'PASS',  # could be unspecified, if so, no pre-filter
                              parseAnno = c('VEP','SNPEFF'))
  # This assumes that the vcf file contains both VEP and SNPEFF output; if not, edit the parseAnno argument.

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
  # check if order somehow messed up
  if(all(vcfs_imp$SAMPLEID == samps_imp$SAMPLEID) & all(rownames(vcfs_imp) == rownames(samps_imp))){
    vcfs_imp %<>% cbind(., samps_imp) #<--- only works if patient / column names are compatible and proper names
  } else {
    warning('Sample order somehow messed up? Not cbinding vcfs with samples.')
  }
  
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
  conCat_veps <- do.call(c, vcfs_imp_l$allVEPmats)
  conCat_snpEffs <- do.call(c, vcfs_imp_l$allSnpEffmats)
  cat('\n Melting list of matrices to long format \n')
  ### setdiff( vcfs_imp_l$vepColnames)
  long_vep_effects_mat <- melt_VEP(veps = conCat_veps, vcf = vcfs_imp, onlyGenes = T, outCols = c('Allele','Consequence', 'IMPACT','SYMBOL','Gene'))  # outCols = common_cols to output all shared VEP fields
  long_snpEff_effects_mat <- melt_VEP(veps = conCat_snpEffs, vcf = vcfs_imp, onlyGenes = T, outCols = c('Allele','Consequence', 'IMPACT','SYMBOL','Gene'))  # outCols = common_cols to output all shared VEP fields
  return(list(vcfs_imp = vcfs_imp, fmats = fmats, long_effects_mat = long_effects_mat, VEPcolnames = common_cols))
}

```


