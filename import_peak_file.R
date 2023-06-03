library(rtracklayer)

# for more info: https://charlesjb.github.io/How_to_import_narrowPeak/
# The rtracklayer package offers ways to import various genomic formats..
#.. such as BED, WIG or GFF/GTF.

# import narrowPeak files
narrow_files <- list.files(pattern = "narrowPeak")

extraCols_narrowPeak <- c(signalValue = "numeric", pValue = "numeric",
                          qValue = "numeric", peak = "integer")

narrowPeak1 <- import(narrow_files[[1]], format = "BED",
                        extraCols = extraCols_narrowPeak)
                        
narrowPeak1

# import broadPeak files
broad_files <- list.files(pattern = "broadPeak")

extraCols_broadPeak <- c(signalValue = "numeric", pValue = "numeric",
                         qValue = "numeric")
                         
broadPeak1 <- import(broad_files[[1]], format = "BED",
                           extraCols = extraCols_broadPeak)

broadPeak1
