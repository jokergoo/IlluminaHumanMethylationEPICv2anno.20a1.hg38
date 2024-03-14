
# another way to get rid of lazy loading is to import the object in into the namespace when
# loading this package
# .onLoad = function(libname, pkgname) {

#     all_vars =  c("IlluminaHumanMethylationEPICv2anno.20a1.hg38",
#                 "Islands.UCSC",
#                 "Locations",
#                 "Manifest",
#                 # "Other",
#                 # "SNPs.Illumina",
#                 "SNPs.141CommonSingle",
#                 "SNPs.142CommonSingle",
#                 "SNPs.144CommonSingle",
#                 "SNPs.146CommonSingle",
#                 "SNPs.147CommonSingle",
#                 "SNPs.150CommonSingle",
#                 "SNPs.151CommonSingle"
#                 )

#     data(list = all_vars, package = pkgname, lib.loc = libname)
#     ns = asNamespace(pkgname)

#     for(var in all_vars) {
#         assign(var, get(var), envir = ns)
#         namespaceExport(ns, var)
#     }
# }


# == title
# Aggregate to the probe-level
#
# == param
# -x A data frame or a matrix.
#
# == details
# If ``x`` is a data frame, it should be one of the following objects: Islands.UCSC, Locations, Other, SNPs.Illumina, SNPs.141CommonSingle, SNPs.142CommonSingle, SNPs.144CommonSingle, 
#  SNPs.146CommonSingle, SNPs.147CommonSingle, SNPs.150CommonSingle, SNPs.151CommonSingle.
#
# If ``x`` is a matrix, it should come from the analysis with the **minfi** package. If multiple rows correspond to the same
# probe, the average value of rows is simply used.
aggregate_to_probes = function(x) {
	nm = as.character(substitute(x))

	if(grepl("(Islands.UCSC|Locations|Other|SNPs.Illumina|SNPs.141CommonSingle|SNPs.142CommonSingle|SNPs.144CommonSingle|SNPs.146CommonSingle|SNPs.147CommonSingle|SNPs.150CommonSingle|SNPs.151CommonSingle)$", nm)) {
		probe_ID = gsub("_.*$", "", rownames(x))
		x = x[!duplicated(probe_ID), , drop = FALSE]
		rownames(x) = gsub("_.*$", "", rownames(x))
		x
	} else if(grepl("Manifest$", nm)) {
		stop("You cannot aggregate to the probes with the `Manifest` object.")
	} else {
		if(all(apply(x, 2, is.numeric))) {
			x = as.matrix(x)
			rn = unique(gsub("_.*$", "", rownames(x)))
			x2 = do.call(rbind, tapply(1:nrow(x), gsub("_.*$", "", rownames(x)), function(ind) {
			    colMeans(x[ind, , drop = FALSE], na.rm = TRUE)
			}, simplify = FALSE))
			x2[rn, , drop = FALSE]
		} else {
			if(!is.null(rownames(x))) {
				probe_ID = gsub("_.*$", "", rownames(x))
				x = x[!duplicated(probe_ID), , drop = FALSE]
				rownames(x) = gsub("_.*$", "", rownames(x))
				x
			} else {
				stop("Wrong format")
			}
		}
	}
}
