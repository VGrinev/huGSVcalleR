#' Build a hash table for the reference genome
#' @description This is a wrapper function that uses the low-level function
#'     buildindex() from package Subread to build the hash table for the
#'     reference genome.
#' @param ref_genome character string giving the name of folder for storing
#'     of index for reference genome.
#' @param ref_fasta character string giving the name of folder for storing
#'     of reference sequences.
#' @param index character string giving the basename of created index files.
#' @param fa charater string giving the name of the FASTA file containing all
#'     the reference sequences.
#' @param memory a numeric value specifying the amount of memory (in megabytes)
#'     used for storing the index during read mapping. Default value is 8000 MB.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return no value is produced but index files are written to the respective
#'     directory.
#' @author Vasily V. Grinev
#' @examples
#' hash <- buildIndexSubread(ref_genome="Reference_Genomes",
#'                           ref_fasta="Files_FASTA",
#'                           index="GRCh38",
#'                           fa="GRCh38.fa",
#'                           memory=8000,
#'                           workDir="/mnt/data/grinev")
#' @export
#' @importFrom Rsubread buildindex
buildIndexSubread <- function(ref_genome,
                              ref_fasta,
                              index,
                              fa,
                              memory=8000,
                              workDir=NULL){
    ### Setting the working directory.
    if (is.null(x=workDir)){
        workDir <- getwd()
    }
    ### Setting the directories.
    refGenome <- paste(workDir, ref_genome, sep="/")
    refFASTA <- paste(workDir, ref_fasta, sep="/")
    ### Calculation of hash table.
    idx <- buildindex(basename=paste(refGenome, index, sep="/"),
                      reference=paste(refFASTA, fa, sep="/"),
                      gappedIndex=TRUE,
                      indexSplit=TRUE,
                      memory=memory,
                      TH_subread=24,
                      colorspace=FALSE)
}
