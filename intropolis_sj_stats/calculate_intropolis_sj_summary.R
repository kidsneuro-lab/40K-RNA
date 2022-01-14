library(data.table)
library(futile.logger)
library(RSQLite)
library(stringr)

flog.appender(appender.file("logs/calculate_intropolis_sj_summary.runtime.log"))
flog.threshold("DEBUG")
flog.info("Starting script")

args <- commandArgs(trailingOnly=TRUE)

if (length(args) > 0) {
  intropolis_sj_file <- args[1]
  intropolis_sj_db_file <- args[2]
  chunk_size <- as.integer(args[3])
} else {
  intropolis_sj_file <- "data/source/intropolis.v1.hg19.tsv"
  intropolis_sj_db_file <- "logs/intropolis_v7_SJ_summary.sqlite"
  chunk_size <- 10000
}


unlink(intropolis_sj_db_file)
intropolis_sj_db <- dbConnect(RSQLite::SQLite(), intropolis_sj_db_file)

flog.info(paste("intropolis SJ file:", intropolis_sj_file))

#### Load intropolis SJ data #### 

flog.info("Load intropolis data")
intropolis_sj_raw <- fread(intropolis_sj_file, header = T, skip = 2, 
                           col.names = c('chr', 'start', 'end', 'strand', 'donor_motif', 'acceptor_motif', 'sample_ids', 'read_counts'),
                           colClasses = c('character',rep('numeric', 2),rep('character', 5)), nThread = 4)


if (nrow(intropolis_sj_raw) > chunk_size) {
  chunks <- trunc(nrow(intropolis_sj_raw) / chunk_size)
  chunks_offsets <- mapply(list, start = c(1, seq(1:chunks) * chunk_size + 1), end = c(seq(1:chunks) * chunk_size, nrow(intropolis_sj_raw)), SIMPLIFY = F)
} else {
  chunks = 1
  chunks_offsets = list(list(start = 1, end = nrow(intropolis_sj_raw)))
}
flog.info(paste("No. of rows:", nrow(intropolis_sj_raw), ", chunks:", chunks))


for (chunk in chunks_offsets) {
  flog.info(paste("--- Processing chunk offset start:", chunk$start, "end:", chunk$end))
  intropolis_sj <- intropolis_sj_raw[chunk$start:chunk$end]
  
  flog.info("Extracting junction locations and assigning sj_key")
  
  intropolis_sj$sj_key = chunk$start:chunk$end

  #### Calculate summary stats across all samples #### 
  flog.info("Calculate summary stats across all samples")
  intropolis_sj_readCounts <- lapply(strsplit(intropolis_sj$read_counts, split = ','), as.numeric)
  
  intropolis_sj$sample_count <- lengths(intropolis_sj_readCounts)
  intropolis_sj$max_sj <- unlist(lapply(intropolis_sj_readCounts, max))
  intropolis_sj$min_sj <- unlist(lapply(intropolis_sj_readCounts, min))
  intropolis_sj$mean_sj <- unlist(lapply(intropolis_sj_readCounts, mean))
  intropolis_sj$median_sj <- unlist(lapply(intropolis_sj_readCounts, median))
  intropolis_sj$q01 <- unlist(lapply(intropolis_sj_readCounts, function(x) { trunc(quantile(x, 0.01)) }))
  intropolis_sj$q05 <- unlist(lapply(intropolis_sj_readCounts, function(x) { trunc(quantile(x, 0.05)) }))
  intropolis_sj$q95 <- unlist(lapply(intropolis_sj_readCounts, function(x) { trunc(quantile(x, 0.95)) }))
  intropolis_sj$q99 <- unlist(lapply(intropolis_sj_readCounts, function(x) { trunc(quantile(x, 0.99)) }))
  
  flog.info("Writing to DB")
  dbWriteTable(intropolis_sj_db, "all_samples_sj_stats", intropolis_sj[,.(chr, start, end, strand, donor_motif, acceptor_motif,
                                                              sample_count, max_sj, min_sj,
                                                              mean_sj, median_sj, q01, q05, q95, q99)], overwrite = FALSE, append = TRUE)
  flog.info("Completed writing to DB")
  
}


dbDisconnect(intropolis_sj_db)














