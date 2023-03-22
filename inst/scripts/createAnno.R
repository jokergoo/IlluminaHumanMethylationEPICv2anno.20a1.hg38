setwd("~/project/development/IlluminaHumanMethylationEPICv2anno.20a1.hg38")


if(basename(getwd()) != "IlluminaHumanMethylationEPICv2anno.20a1.hg38") {
     stop("Working dir should be 'IlluminaHumanMethylationEPICv2anno.20a1.hg38' (in the package dir).")
}

library(minfi)

manifestFile <- "../data_EPIC_v2/EPIC-8v2-0_A1.csv"

source("inst/scripts/read.manifest.EPIC.R")
maniTmp <- read.manifest.EPIC(manifestFile)
anno <- maniTmp$manifest

## Checking start
library(illuminaio)
epic <- readIDAT("../data_EPIC_v2/206891110001_R01C01_Grn.idat")
address.epic <- as.character(epic$MidBlock)
dropCpGs <- anno$Name[anno$AddressB != "" & !anno$AddressB %in% address.epic]
dropCpGs <- anno$Name[anno$AddressA != "" & !anno$AddressA %in% address.epic]
table(substr(dropCpGs, 1,2))
## Checking done

# probe names can be duplicated and only affect the following three columns
# IlmnID = tapply(anno$IlmnID, anno$Name, paste, collapse = ";")
# AddressA_ID = tapply(anno$AddressA_ID, anno$Name, paste, collapse = ";")
# AddressB_ID = tapply(anno$AddressB_ID, anno$Name, paste, collapse = ";")

# anno = anno[!duplicated(anno$Name), ]
# anno$IlmnID = IlmnID[anno$Name]
# anno$AddressA_ID = AddressA_ID[anno$Name]
# anno$AddressB_ID = AddressB_ID[anno$Name]

anno = anno[anno$CHR %in% paste0("chr", c(1:22, "X", "Y")), ]
manifestList <- maniTmp$manifestList


## Manifest package
IlluminaHumanMethylationEPICv2manifest <- do.call(IlluminaMethylationManifest,
                                                list(TypeI = manifestList$TypeI,
                                                     TypeII = manifestList$TypeII,
                                                     TypeControl = manifestList$TypeControl,
                                                     TypeSnpI = manifestList$TypeSnpI,
                                                     TypeSnpII = manifestList$TypeSnpII,
                                                     annotation = "IlluminaHumanMethylationEPICv2"))
## Annotation package
anno$IlmnID <- NULL
nam <- names(anno)
names(nam) <- nam
nam[c("AddressA_ID", "AddressB_ID", "AlleleA_ProbeSeq", "AlleleB_ProbeSeq",
            "Infinium_Design_Type", "Next_Base", "Color_Channel")] <-  c("AddressA", "AddressB",
                                                                         "ProbeSeqA", "ProbeSeqB",
                                                                         "Type", "NextBase", "Color")

names(nam) <- NULL
names(anno) <- nam
rownames(anno) <- anno$Name
anno <- anno[intersect(anno$Name, unique(getManifestInfo(IlluminaHumanMethylationEPICv2manifest, type = "locusNames"))),]

Locations <- anno[, c("CHR", "MAPINFO")]
names(Locations) <- c("chr", "pos")
Locations$pos <- as.integer(Locations$pos)
Locations$strand <- ifelse(anno$Strand_FR == "F", "+", "-")
table(Locations$chr, exclude = NULL)
rownames(Locations) <- anno$Name
Locations <- as(Locations, "DataFrame")

Manifest <- anno[, c("Name", "AddressA", "AddressB",
                     "ProbeSeqA", "ProbeSeqB", "Type", "NextBase", "Color")]
Manifest <- as(Manifest, "DataFrame")

Islands.UCSC <- anno[, c("UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island")]
names(Islands.UCSC) <- c("Islands_Name", "Relation_to_Island")
Islands.UCSC <- as(Islands.UCSC, "DataFrame")
Islands.UCSC$Relation_to_Island[Islands.UCSC$Relation_to_Island == ""] <- "OpenSea"
Islands.UCSC$Relation_to_Island[grepl("Island", Islands.UCSC$Relation_to_Island)] = "Island"
Islands.UCSC$Relation_to_Island[grepl("Shore", Islands.UCSC$Relation_to_Island)] = "Shore"
Islands.UCSC$Relation_to_Island[grepl("Shelf", Islands.UCSC$Relation_to_Island)] = "Shelf"
table(Islands.UCSC$Relation_to_Island, exclude = NULL)

SNPs.Illumina <- anno[, c("SNP_ID", "SNP_DISTANCE", "SNP_MinorAlleleFrequency")]
SNPs.Illumina <- as(SNPs.Illumina, "DataFrame")

usedColumns <- c(names(Manifest), names(SNPs.Illumina), 
                 c("CHR", "MAPINFO", "Strand",
                   "Chromosome_36", "Coordinate_36", "Genome_Build"),
                 c("UCSC_CpG_Islands_Name", "Relation_to_UCSC_CpG_Island"))
Other <- anno[, setdiff(names(anno), c(usedColumns, "Source_Seq", "Forward_Sequence", "Top_Sequence", "GencodeV41_Accession",
                                       "UCSC_RefGene_Group", "UCSC_RefGene_Accession", "GencodeV41_Group", "Regulatory_Feature_Name",
                                       "DNase_Hypersensitivity_NAME", "Encode_CisReg_Site", "Encode_CisReg_Site_Evid",
                                       "OpenChromatin_NAME", "OpenChromatin_Evidence_Count"))]
nam <- names(Other)
nam <- sub("_NAME", "_Name", nam)
nam[nam == "X450k_Enhancer"] <- "Methyl450_Enhancer"
nam
names(Other) = nam
Other <- as(Other, "DataFrame")

## We now use an exisitng grSnp object containing a GRanges of relevant SNPs.
## This is created in a separate script

##
## SNP overlap
##

map <- cbind(Locations, Manifest)
map <- GRanges(seqnames = map$chr, ranges = IRanges(start = map$pos, width = 1),
               Strand = map$strand, Type = map$Type)
map <- minfi:::.getProbePositionsDetailed(map)
names(map) <- rownames(Locations)

## dbSNP

system("
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp141Common.txt.gz
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp142Common.txt.gz
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp144Common.txt.gz
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp146Common.txt.gz
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp147Common.txt.gz
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp150Common.txt.gz
curl -O https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/snp151Common.txt.gz

gunzip -c snp141Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp141Common_small.txt.gz
gunzip -c snp142Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp142Common_small.txt.gz
gunzip -c snp144Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp144Common_small.txt.gz
gunzip -c snp146Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp146Common_small.txt.gz
gunzip -c snp147Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp147Common_small.txt.gz
gunzip -c snp150Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp150Common_small.txt.gz
gunzip -c snp151Common.txt.gz | cut -f2,3,4,5,7,8,12,25 | gzip - > snp151Common_small.txt.gz
")

file.remove("snp141Common.txt.gz")
file.remove("snp142Common.txt.gz")
file.remove("snp144Common.txt.gz")
file.remove("snp146Common.txt.gz")
file.remove("snp147Common.txt.gz")
file.remove("snp150Common.txt.gz")
file.remove("snp151Common.txt.gz")

processUCSCsnp <- function(snpfile) {
    require(GenomicRanges)
    cat("Reading file\n")
    df <- read.delim(gzfile(snpfile), header = FALSE,
                     stringsAsFactors = FALSE)
    names(df) <- c("chr", "start", "end", "name", "strand",
                   "refNCBI", "class", "alleleFreqs")
    print(table(df$chr))
    cat("Only keeping chrs 1-22, X, Y\n")
    df <- df[df$chr %in% paste0("chr", c(1:22, "X", "Y")),]
    print(table(df$class))
    cat("Only keeping class 'single'\n")
    df <- df[df$class == "single",]
    cat("Computing MAF\n")
    df$alleleFreqs <- sub(",$", "", df$alleleFreqs)
    sp <- strsplit(df$alleleFreqs, ",")
    minFreq <- sapply(sp, function(xx) min(as.numeric(xx)))
    cat("Instantiating object\n")
    grSNP <- GRanges(seqnames = df$chr, strand = df$strand,
                     ranges = IRanges(start = df$start + 1, end = df$end),
                     MAF = minFreq, ref = df$refNCBI)
    names(grSNP) <- df$name
    grSNP
}


grSnp141CommonSingle <- processUCSCsnp("snp141Common_small.txt.gz")
grSnp142CommonSingle <- processUCSCsnp("snp142Common_small.txt.gz")
grSnp144CommonSingle <- processUCSCsnp("snp144Common_small.txt.gz")
grSnp146CommonSingle <- processUCSCsnp("snp146Common_small.txt.gz")
grSnp147CommonSingle <- processUCSCsnp("snp147Common_small.txt.gz")
grSnp150CommonSingle <- processUCSCsnp("snp150Common_small.txt.gz")
grSnp151CommonSingle <- processUCSCsnp("snp151Common_small.txt.gz")

file.remove("snp141Common_small.txt.gz")
file.remove("snp142Common_small.txt.gz")
file.remove("snp144Common_small.txt.gz")
file.remove("snp146Common_small.txt.gz")
file.remove("snp147Common_small.txt.gz")
file.remove("snp150Common_small.txt.gz")
file.remove("snp151Common_small.txt.gz")

SNPs.141CommonSingle <- minfi:::.doSnpOverlap(map, grSnp141CommonSingle)
SNPs.142CommonSingle <- minfi:::.doSnpOverlap(map, grSnp142CommonSingle)
SNPs.144CommonSingle <- minfi:::.doSnpOverlap(map, grSnp144CommonSingle)
SNPs.146CommonSingle <- minfi:::.doSnpOverlap(map, grSnp146CommonSingle)
SNPs.147CommonSingle <- minfi:::.doSnpOverlap(map, grSnp147CommonSingle)
SNPs.150CommonSingle <- minfi:::.doSnpOverlap(map, grSnp150CommonSingle)
SNPs.151CommonSingle <- minfi:::.doSnpOverlap(map, grSnp151CommonSingle)

annoNames <- c("Locations", "Manifest", "SNPs.Illumina",
               "SNPs.141CommonSingle", "SNPs.142CommonSingle", "SNPs.144CommonSingle",
               "SNPs.146CommonSingle", "SNPs.147CommonSingle", "SNPs.150CommonSingle",
               "SNPs.151CommonSingle", "Islands.UCSC", "Other")
for(nam in annoNames) {
    cat(nam, "\n")
    save(list = nam, file = file.path("data", paste(nam, "rda", sep = ".")), compress = "xz")
}

annoStr <- c(array = "IlluminaHumanMethylationEPICv2",
             annotation = "20a1",
             genomeBuild = "hg38")
defaults <- c("Locations", "Manifest", "SNPs.151CommonSingle", "Islands.UCSC", "Other")
pkgName <- sprintf("%sanno.%s.%s", annoStr["array"], annoStr["annotation"],
                    annoStr["genomeBuild"])

annoObj <- IlluminaMethylationAnnotation(objectNames = annoNames, annotation = annoStr,
                              defaults = defaults, packageName = pkgName)

assign(pkgName, annoObj)
save(list = pkgName,
     file = file.path("data", paste(pkgName, "rda", sep = ".")), compress = "xz")
# sessionInfo()




