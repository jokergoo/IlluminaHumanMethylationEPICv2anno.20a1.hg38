
library(Biostrings)
read.manifest.EPIC <- function(file) {
    # NOTE: As is, requires grep
    control.line <- system(
        sprintf("grep -n \\\\[Controls\\\\] %s", file), intern = TRUE)
    control.line <- as.integer(sub(":.*", "", control.line))
    stopifnot(length(control.line) == 1 &&
                  is.integer(control.line) &&
                  !is.na(control.line))
    assay.line <- system(
        sprintf("grep -n \\\\[Assay\\\\] %s", file), intern = TRUE)
    assay.line <- as.integer(sub(":.*", "", assay.line))
    stopifnot(length(assay.line) == 1 &&
                  is.integer(assay.line) &&
                  !is.na(assay.line))
    colNames <- readLines(file, n = assay.line + 1L)[assay.line + 1L]
    colNames <- strsplit(colNames, ",")[[1]]
    colClasses <- rep("character", length(colNames))
    names(colClasses) <- colNames
    names(colClasses) <- make.names(names(colClasses))
    colClasses[c("MAPINFO")] <- "integer"
    manifest <- read.table(
        file = file,
        header = FALSE,
        col.names = names(colClasses),
        sep = ",",
        comment.char = "",
        quote = "",
        skip = assay.line + 1L,
        colClasses = colClasses,
        nrows = control.line - assay.line - 2L)
    manifest$Name = manifest$IlmnID
    manifest$AddressA_ID <- gsub("^0*", "", manifest$AddressA_ID)
    manifest$AddressB_ID <- gsub("^0*", "", manifest$AddressB_ID)
    TypeI <- manifest[
        manifest$Infinium_Design_Type == "I",
        c("Name", "AddressA_ID", "AddressB_ID", "Color_Channel", "Next_Base",
          "AlleleA_ProbeSeq", "AlleleB_ProbeSeq")]
    names(TypeI)[c(2, 3, 4, 5, 6, 7)] <- c(
        "AddressA", "AddressB", "Color", "NextBase", "ProbeSeqA", "ProbeSeqB")
    TypeI <- as(TypeI, "DataFrame")
    TypeI$ProbeSeqA <- DNAStringSet(TypeI$ProbeSeqA)
    TypeI$ProbeSeqB <- DNAStringSet(TypeI$ProbeSeqB)
    TypeI$NextBase <- DNAStringSet(TypeI$NextBase)
    TypeI$nCpG <- as.integer(
        oligonucleotideFrequency(TypeI$ProbeSeqB, width = 2)[, "CG"] - 1L)
    TypeI$nCpG[TypeI$nCpG < 0] <- 0L
    TypeSnpI <- TypeI[grep("^rs", TypeI$Name), ]
    TypeI <- TypeI[-grep("^rs", TypeI$Name), ]

    TypeII <- manifest[
        manifest$Infinium_Design_Type == "II",
        c("Name", "AddressA_ID", "AlleleA_ProbeSeq")]
    names(TypeII)[c(2,3)] <- c("AddressA", "ProbeSeqA")
    TypeII <- as(TypeII, "DataFrame")
    TypeII$ProbeSeqA <- DNAStringSet(TypeII$ProbeSeqA)
    TypeII$nCpG <- as.integer(letterFrequency(TypeII$ProbeSeqA, letters = "R"))
    TypeII$nCpG[TypeII$nCpG < 0] <- 0L
    TypeSnpII <- TypeII[grep("^rs", TypeII$Name), ]
    TypeII <- TypeII[-grep("^rs", TypeII$Name), ]

    controls <- read.table(
        file = file,
        skip = control.line,
        sep = ",",
        comment.char = "",
        quote = "",
        colClasses = c(rep("character", 5)))[, 1:5]
    TypeControl <- controls[, 1:4]
    names(TypeControl) <- c("Address", "Type", "Color", "ExtendedType")
    TypeControl <- as(TypeControl, "DataFrame")

    list(
        manifestList = list(
            TypeI = TypeI,
            TypeII = TypeII,
            TypeControl = TypeControl,
            TypeSnpI = TypeSnpI,
            TypeSnpII = TypeSnpII),
        manifest = manifest,
        controls = controls)
}
