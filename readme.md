IlluminaHumanMethylationEPICv2anno.20a1.hg38
-----------------------------------

This package provides annotations for Illumina methylation EPIC array v2.0. The data is based on the file https://support.illumina.com/content/dam/illumina-support/documents/downloads/productfiles/methylationepic/MethylationEPIC%20v2%20Files.zip from 
  https://support.illumina.com/array/array_kits/infinium-methylationepic-beadchip-kit/downloads.html.

When using with the **minfi** package, you can manually set the "annotation" element by providing the suffix (removing the string "IlluminaHumanMethylationEPICv2anno."):

```r
RGset = read.metharray.exp(...)

# explained in the IlluminaHumanMethylationEPICv2manifest package
annotation(RGset)["array"] = "IlluminaHumanMethylationEPICv2"

annotation(RGset)["annotation"] = "20a1.hg38"
```


## Install

```r
BiocManager::install("jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38")
# or
remotes::install_github("jokergoo/IlluminaHumanMethylationEPICv2anno.20a1.hg38")
```


## Licence

Artistic-2.0

