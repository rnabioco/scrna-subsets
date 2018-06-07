
# compression scheme used by 10x, function below is an internal function taken from their package
#' Compress a vector of sequences into 2-bit representation
#'
#' Compress a vector of sequences into integer64 objects containing
#' a 2-bit representation. Ns are not allowed
#' @param seqs Vector of nucleotide sequences to compress
#' @return A vector of integer64 objects
compress_sequences <- function(seqs) {
  if (any(grepl('[^ACGT]', seqs))) {
    stop("At least one sequence contains Ns")
  }
  nuc_to_int <- as.integer(0:3)
  names(nuc_to_int) <- c('A', 'C', 'G', 'T')
  
  chars <- do.call(rbind, strsplit(seqs, ''))
  nuc_ints <- matrix(nuc_to_int[chars], nrow=length(seqs))
  result <- integer64(length(seqs))
  for(i in 1:ncol(nuc_ints)) {
    result <- result * as.integer(4) + nuc_ints[,i]
  }
  
  result
}

#write a custom decompressing function in R

#' decompress the 10x encoded barcodes. 
#' @param seqs_int64 Vector of int64 2bit encoded sequences to decompress
#' @param lens Length of compressed barcode (int)
#' @return A vector of DNA sequences
#' @example 
#' seqs <- c("ATCG", "GGAT", "CTGA", "TGAC")
#' seqs_binary <- compress_sequences(seqs)
#' out_seqs <- decompress_sequences(seqs_binary, 4)
#' all(seqs == out_seqs)
decompress_sequences <- function(seqs, lens){
  # hacky and takes about 15 seconds per 1e6  sequences
  bit_code <- c('A', 'C', 'G', 'T')
  names(bit_code) <- c("00", "01", "10", "11")
  
  b <- as.bitstring(seqs)
  b <- do.call(rbind, strsplit(b, ""))
  nbit <- lens * 2
  nbit_pos <- 64 - (nbit - 1)
  b <- b[, nbit_pos:64]
  pos <- seq(1, nbit, by = 2)
  result = ""
  for(i in pos){
    x <- paste0(b[, i], b[, (i + 1)])
    x_n <- bit_code[x]
    result <- paste0(result, x_n)
  }
  result
}



