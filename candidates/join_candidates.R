samples = c('cortex_bulk', 'heart_bulk', 'neuron_2', 'neuron_3', 'neuron_51',
            'neuron_6', 'neurons_100batch')

# Columns that must exist in the file, used for joining
join_cols = c("chr", "start", "end", "ref_len", "unit", "region")
data_cols = c("raw_alleles", "call", "genotype", "pval", "allele_summaries")

reader <- function(file) {
    cat(sprintf("Reading %s...\n", file))
    read.table(file, sep="\t", header=TRUE)
}

# Make an empty data.frame with the expected columns
df <- do.call(data.frame, as.list(c(join_cols, data_cols)))  # Dummy row
colnames(df) <- c(join_cols, data_cols)
df <- df[FALSE,]  # Remove the row of dummy values


all <- TRUE
for (sample in samples) {
    t <- reader(sprintf('%s/%s.flank10.minmapq36_minunits4_minsupp15_maxrefdiff80.candidates.txt', sample, sample))

    df <- merge(df, t, sort=FALSE, by=join_cols, all=all,
                suffixes=c('', paste('.', sample, sep='')))
    all <- FALSE
}

# All of the data columns without a .XXX suffix are the original dummy
# data columns, which are all NA.
for (col in data_cols) {
    df[col] <- NULL
}

write.table(df, file="all_called.flank10.minmapq36_minunits4_minsupp15_maxrefdiff80.candidates.txt",
            quote=FALSE, sep="\t", row.names=FALSE)
