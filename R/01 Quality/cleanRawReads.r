## FASTQ raw data preprocessing function.
preproFASTQ = function(dir_fastq, fastqFile){
#  Loading of the FASTQ file.
raw.data = readFastq(dirPath = dir_fastq, pattern = fastqFile, withIds = TRUE)
#  Data preprocessing, step 1 (trimming of the low quality ends).
raw.data = narrow(x = raw.data, start = 21, end = 120, width = NA, use.names = TRUE)
cat("Low quality ends of reads have been trimmed.\n")
flush.console()
#  Data preprocessing, step 2 (removing of the reads with ambiguous nucleotides).
raw.data = clean(raw.data)
cat("Reads with ambiguous nucleotides have been removed from the FASTQ file.\n")
flush.console()
#  Data preprocessing, step 3 (removing of the low quality reads).
raw.data = raw.data[rowMeans(as(quality(raw.data), "matrix")) > 29]
cat("Low quality reads have been removed from the FASTQ file.\n")
flush.console()
#  Data preprocessing, step 4 (removing of the duplicated reads).
raw.data = raw.data[occurrenceFilter(min = 1L,
                                     max = 1L,
                                     withSread = c(TRUE),
                                     duplicates = c("sample"))(raw.data)]
cat("Duplicated reads have been removed from the FASTQ file.\n")
flush.console()
return(raw.data)
}
