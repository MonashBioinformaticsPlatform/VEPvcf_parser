**packages required:**

mutSignatures

stringr

magrittr

----

**Useage example:**

get all .vcf.gz files in working directory and parse into a table (filtering for FILTER='PASS' for speed).

```
library(stringr)
library(mutSignatures)
library(magrittr)
source('VEPvcf_parser.R')

vcfFNs <- dir(pattern = '.+\\.vcf\\.gz$')
sampnames <- vcfFNs %>% str_remove(pattern='VEP.ann.vcf.gz')  # assuming vcf files are named this way
vcfs <- importVEPVCFfiles(vcfFiles =  paste0(indir, '/', vcfFNs), sampleNames = sampnames, filterIn = 'PASS')
```
The output is a dataframe, is one line per variant, with a matrix of variant effects in each cell of column `VEP_matrix`

----

To convert list of matrices in column `VEP_matrix` to long format (ignoring effects that are not explicitly linked to a gene) use `melt_VEP`:
```
vcfs$HANDLE <- paste0(vcfs$CHROM,'_', vcfs$POS, '_', vcfs$REF, '_', vcfs$ALT) # make a unique handle string for each variant
long <- melt_VEP(vcfs,handle = 'HANDLE' , outCols=c('Allele','Consequence', 'IMPACT','SYMBOL','Gene'))
```
Output is a dataframe, one line per **variant / affected gene** combination, with the columns from the original that are nominated in the `outCols` argument.
