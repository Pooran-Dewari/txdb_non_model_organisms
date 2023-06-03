#START : part1: txdb tutorial ..............................................
# txdb objects are containers to store transcript annotations

library(GenomicFeatures) # makeTxDbFromGFF; loadDb

#create a TxDb object from gff file, need to do it only once
TxDb_Ss <- makeTxDbFromGFF("Salmo_salar.Ssal_v3.1.106.gff3",
                           dataSource = "Ensembl", organism = "Salmo salar")

#save TxDb for future use
saveDb(TxDb_Ss, file = "TxDb_Salmo_salar")

# for now, just some basic playing around txdb
txdb <- TxDb_Ss
txdb # info about the txdb object
columns(txdb)
keytypes(txdb)
glimpse_txdb <- head(select(txdb, columns=columns(txdb), keys=keys(txdb), keytype=c("GENEID"))) #glimpse into everything in txdb
glimpse_txdb #can see that GENEID is ensembl gene ID

# we can use select method to retrieve different columns for custom geneid query
my_genes <- c("ENSSSAG00000000040", "ENSSSAG00000000057") # keys are GENEID
cols <- c("CDSCHROM", "CDSSTART", "CDSEND")
my_genes_df <- select(txdb, keys = my_genes, keytype="GENEID", columns=cols) #retrieve columns and save to df

#returning GRanges objects
#the most common operations for a TxDb object is to retrieve the
#...genomic coordinates or ranges for exons, transcripts or coding sequences.
#let's get promoters coordinates
PR_salmo_salar <- promoters(txdb, upstream=2000, downstream=500)

# return all exons for the organism
EX_salmo_salar <- exons(txdb)
#cds(txdb)
#transcripts(txdb)
#genes(txdb)


# load txdb and select seqlevels, ie. chromosome
txdb <- loadDb(file = "TxDb_Salmo_salar")
seqlevels(txdb) #this will show chr for txdb
seqlevels(txdb) <- c("1","2","3", "4", "5", "6", "7", "8", "9", "10",
                     "11", "12", "13", "14", "15", "16", "17", "18", "19",
                     "20", "21", "22", "23", "24", "25", "26", "27", "28", "29") #restrict to chr1 to 29 
#seqlevels(txdb) <- seqlevels0(txdb)     #use this to reset to default chr

#END : part1: txdb tutorial ..............................................
