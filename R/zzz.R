
# another way to get rid of lazy loading is to import the object in into the namespace when
# loading this package
# .onLoad = function(libname, pkgname) {

#     all_vars =  c("IlluminaHumanMethylationEPICv2anno.20a1.hg38",
#                 "Islands.UCSC",
#                 "Locations",
#                 "Manifest",
#                 "Other",
#                 "SNPs.Illumina",
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
