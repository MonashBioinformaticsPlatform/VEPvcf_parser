**packages required:**

mutSignatures

stringr

magrittr

----

**example:**

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



