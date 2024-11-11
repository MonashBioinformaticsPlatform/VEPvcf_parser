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

The output after melting is a dataframe, one line per **variant / affected gene / effect** combination, with the columns from the original that are nominated in the `outCols` argument. Therefore, the number of rows will be greater than in the original VCF file. To refer back to the original entry, the new dataframe will contain an `original_rowname` column which refers back to the rownames of the vcf in the list `vcf_in` (output of the vcf2list() function) or, if that is not provided, the names of the elements in vepMatrices. This behaviour can be over-ridden by providing original rownames as a vector of unique names (same length as `veps`) using the `handles` argument.  

### Limitations:

This is only tested on one-alt-allele-per-line vcfs, not comma-separated alt alleles.
If you have multiallelic entries on single rows, use "bcftools norm -m - " (note the minus sign) to split multiallelic sites over several rows.
Also it expects some format fields that might not be present in some vcfs ('DP', 'GQ', 'GT' and 'AD'), if these are absent then you will need to specify `formFields` argument with vaules that can be accepted by `vcfR::extract.gt()`. 

### Do lots of files at once:

If vcf files are from non-overlapping genomic regions, with the same set of samples (e.g. results from chunking a file over genomic regions), can use `readMultiVcfs`.

See the associated script `example_VEPvcf_parser.R` for details.


