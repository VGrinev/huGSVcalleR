# This function is completely identical to the internal function cutRseq
# from FastqCleaner package. It is based on trimLRpattern method from Biostrings.
# 

# Arguments:

# subject - DNAString or DNAStringSet object;
# Rpattern - 3’ pattern, DNAString object;
# with.indels -  Allow indels?;
# fixed - Parameter passed to trimLRPatterns Default ’subject’, 
#         ambiguities in the pattern only are interpreted as wildcard. 
#         See the argument fixed in trimLRPatterns;
# error_rate - Error rate (value in [0, 1]). The error rate is the proportion of
#         mismatches allowed between the adapter and the aligned portion of the 
#         subject. For a given adapter A, the number of allowed mismatches between 
#         each subsequence s of A and the subject is computed as: error_rate * L_s, 
#         where L_s is the length of the subsequence s;
# anchored - Can the adapter or partial adapter be within the sequence? 
#         (anchored = FALSE) or only in the terminal regions of the sequence? 
#         (anchored = TRUE). Default TRUE (trim only flanking regions);
# ranges - return ranges (default FALSE);
# checks - Perform internal checks? Default TRUE;
# min_match_flank - Do not trim in flanks of the subject, if a match has 
#                   min_match_flank of less length. Default 1L 
#                   (only trim with >=2 coincidences in a flank match)

cutRigthAdapter <- function(subject, Rpattern, 
                            with.indels = FALSE,
                            fixed = "subject", 
                            error_rate = 0.2, 
                            anchored = TRUE, 
                            checks = TRUE, 
                            min_match_flank = 2L, 
                            ...) {
  
  Rpattern <- DNAString(Rpattern)
  
  if (error_rate > 1 || error_rate < 0) {
    stop("error_rate must be a number between 0 and 1")
  }
  
  if (checks) {
    
    if (!is(Rpattern, "DNAString")) {
      stop("Rpattern must be a character string or a DNAString object")
    }
    
    csub <- class(subject)
    if (csub != "DNAStringSet") {
      stop("subject must be a DNAString or DNAStringSet object")
    }
    
    if (csub == "DNAString") {
      subject <- as(subject[[1]], "DNAStringSet")
    }
    
  }
  
  p <- length(Rpattern)
  s_width <- width(subject)
  s <- max(width(subject))
  
  if(error_rate > 0) {
    flank_seq <- as.integer(seq_len(p) * error_rate)
  } else {
    flank_seq <- rep(0, length(seq_len(p)))
  }
  
  
  if (min_match_flank >= 1L) {
    if (p > min_match_flank) {
      flank_seq[seq_len(min_match_flank)] <- -1
    } else {
      return(subject)
    }
  }
  
  if(!anchored) {
    Rpattern <- as.character(Rpattern)
    maxlen <-  max(width(subject)) - nchar(Rpattern)
    if(maxlen > 0) {
      Rpattern <- paste0(Rpattern,  paste(rep("N",maxlen), collapse = ""))
    }
    Rpattern <- DNAString(Rpattern)
    flank_seq <- c(flank_seq, rep(0,maxlen))
  }
  
  out <- trimLRPatterns(Rpattern = Rpattern,
                        subject = subject, 
                        max.Rmismatch = flank_seq, 
                        with.Rindels = with.indels, 
                        Rfixed = fixed, 
                        ...)
  if (ranges) {
    out <- IRanges::IRanges(start = rep(1, length(out)), end = width(out))
  }
  
  out
}