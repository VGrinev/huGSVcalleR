#' Call of single nucleotide variations using HaplotypeCaller GATK
#' @description Internal wrapper around GATK HaplotypeCaller for SNV
#'     calling. Called from callSNVs(criteria = "gatk").
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character string giving the name of input BAM file.
#' @param postfix character string that will be added to the output VCF file
#'     name as postfix. NA by default at which the postfix "GATK based SNVs"
#'     will be used.
#' @param fastaDir character string giving the name (or path to and name) of
#'     folder for storing of fai-indexed reference genome.
#' @param faFile charater string giving the name of file containing reference
#'     genome in format FASTA or FA.
#' @param intervals character string specifying the name of the tab-delimited
#'     TXT file containing coordinates of genomic location(-s) of interest. The
#'     default value is NULL. If so, this file must contains the following four
#'     fields:
#'     i) seqnames - the name of chromosome or scaffold with prefix "chr";
#'     ii) start   - start coordinate of genomic interval of interest;
#'     iii) end    - end coordinate of genomic interval of interest;
#'     iv) strand  - (optionally) strand information about genomic interval
#'                   of interest.
#'     TXT file with genomic location(-s) must be in working directory.
#' @param pcr_model GATK argument --pcr-indel-model. Possible values: NONE,
#'     AGRESSIVE, HOSTILE, CONSERVATIVE (by default). See HaplotypeCaller
#'     at https://gatk.broadinstitute.org for further details.
#' @param soft_clip GATK argument --dont-use-soft-clipped. TRUE by default.
#' @param ERC GATK ERC mode. "NONE" by default. Possible values: NONE, GVCF or
#'     BP_RESOLUTION. See HaplotypeCaller at https://gatk.broadinstitute.org
#'     for further details.
#' @param cores number of CPU cores to use for Hidden Markov Models. If "MAX"
#'     (by default), all available cores will be detected and used.
#' @param SE if TRUE (by default), standart error determined by GATK will be
#'     writen in file STDERR.HC.TXT.
#' @param gatk_path character string giving the path to internal or external
#'     (system) GATK executive script. This path can be obtained by running
#'     Find.Gatk() or settled manually. Default value is "gatk", i.e. function
#'     will uses GATK systemly from $PATH.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#'     fastaDir, faFile, bamDir and bamFiles must be somewhere inside workDir.
#' @details HaplotypeCaller performs SNVs and short indels calling by local de
#'     novo assembly of mapped reads. GATK is not natively available on
#'     Windows; on Windows, use WSL2 or Docker and run GATK manually. See the
#'     INSTALL file for details.
#' @return a new VCF file with detected SNVs is written to \code{bamDir}.
#'     The function is called for its side effect and returns \code{NULL}
#'     invisibly.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @keywords internal
#' @noRd
#' @importFrom parallel detectCores
gatkCalleR <- function(bamDir = NULL,
    bamFile,
    postfix = NA,
    fastaDir = NULL,
    faFile,
    intervals = NULL,
    pcr_model = "CONSERVATIVE",
    soft_clip = TRUE,
    ERC = "NONE",
    cores = "MAX",
    SE = TRUE,
    gatk_path = NULL,
    workDir = NULL) {
    ### Platform check.
    if (.Platform$OS.type == "windows") {
        stop(
            "GATK is not natively supported on Windows. ",
            "Please use WSL2 or Docker and run GATK manually, ",
            "then use callSNVs(criteria = 'fisher') or another method. ",
            "See the INSTALL file for details."
        )
    }
    ### Validation of GATK path.
    .validateGatkPath(gatk_path)
    ### Setting the working directory.
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    ### Full paths to input files.
    path1 <- .resolveBamPath(bamDir, workDir)
    path2 <- .resolveFastaPath(fastaDir, workDir)
    ### Building the argument list.
    args <- .buildHaplotypeCallerArgs(
        path1, bamFile, path2, faFile, postfix, intervals,
        pcr_model, soft_clip, ERC, cores, SE, workDir
    )
    ### Running the command.
    .runHaplotypeCaller(gatk_path, args, path1, bamFile, SE)
    ### Returning nothing.
    invisible(NULL)
}

#' Validate GATK executable path
#' @description Internal helper for gatkCalleR(). Checks that the GATK
#'     executable is specified and exists either as a file or in $PATH.
#' @param gatk_path character string with the GATK executable path, or NULL.
#' @return invisible(NULL). Stops on error.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @keywords internal
#' @noRd
.validateGatkPath <- function(gatk_path) {
    if (is.null(x = gatk_path)) {
        stop(
            "Cannot run GATK without path. Please provide full path to GATK ",
            "executive script, or 'gatk' if it can be called from $PATH."
        )
    }
    if (!file.exists(gatk_path) && !nzchar(Sys.which(gatk_path))) {
        stop(
            "GATK executable not found: ", gatk_path,
            ". Provide a valid path or ensure 'gatk' is in $PATH."
        )
    }
    invisible(NULL)
}

#' Resolve BAM directory path
#' @description Internal helper for gatkCalleR(). Returns the full path to
#'     the BAM directory.
#' @param bamDir character string with the BAM directory, or NULL.
#' @param workDir character string with the working directory.
#' @return character string with the full path.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @keywords internal
#' @noRd
.resolveBamPath <- function(bamDir, workDir) {
    if (is.null(x = bamDir)) workDir else paste(workDir, bamDir, sep = "/")
}

#' Resolve FASTA directory path
#' @description Internal helper for gatkCalleR(). Returns the full path to
#'     the FASTA directory.
#' @param fastaDir character string with the FASTA directory, or NULL.
#' @param workDir character string with the working directory.
#' @return character string with the full path.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @keywords internal
#' @noRd
.resolveFastaPath <- function(fastaDir, workDir) {
    if (is.null(x = fastaDir)) workDir else paste(workDir, fastaDir, sep = "/")
}

#' Build HaplotypeCaller argument list
#' @description Internal helper for gatkCalleR(). Builds the character vector
#'     with all GATK HaplotypeCaller arguments.
#' @param path1 character string with the BAM directory.
#' @param bamFile character string with the BAM file name.
#' @param path2 character string with the FASTA directory.
#' @param faFile character string with the FASTA file name.
#' @param postfix character string with the output postfix, or NA.
#' @param intervals character string with the intervals file, or NULL.
#' @param pcr_model character string with the PCR indel model.
#' @param soft_clip logical, use soft-clipped bases?
#' @param ERC character string with the ERC mode.
#' @param cores number or "MAX".
#' @param SE logical, save standard error?
#' @param workDir character string with the working directory.
#' @return character vector with GATK arguments.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildHaplotypeCallerArgs <- function(path1, bamFile, path2, faFile,
    postfix, intervals, pcr_model,
    soft_clip, ERC, cores, SE, workDir) {
    input <- paste("-I", file.path(path1, bamFile), sep = " ")
    bam_basename <- gsub(pattern = ".bam", replacement = "", x = bamFile)
    out_name <- if (!is.na(x = postfix)) {
        paste(bam_basename, postfix, "vcf", sep = ".")
    } else {
        paste(bam_basename, "HCaller", "vcf", sep = ".")
    }
    output <- paste("-O", file.path(path1, out_name), sep = " ")
    fa <- paste("-R", file.path(path2, faFile), sep = " ")
    if (!is.null(x = intervals)) {
        fa <- paste(fa, .buildIntervalsArgument(intervals, workDir))
    }
    pcr <- paste("--pcr-indel-model", pcr_model)
    softclip <- paste("--dont-use-soft-clipped-bases", soft_clip)
    erc <- paste("-ERC", ERC)
    if (cores == "MAX") {
        cores <- detectCores()
    }
    cores <- paste("--native-pair-hmm-threads", cores)
    c("HaplotypeCaller", input, output, fa, pcr, softclip, erc, cores)
}

#' Build intervals argument for GATK
#' @description Internal helper for gatkCalleR(). Builds the --intervals
#'     argument from a tab-delimited TXT file.
#' @param intervals character string with the intervals file name.
#' @param workDir character string with the working directory.
#' @return character string with the --intervals argument.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildIntervalsArgument <- function(intervals, workDir) {
    GRs <- read.table(
        file = paste(workDir, intervals, sep = "/"),
        sep = "\t",
        header = TRUE,
        quote = "\"",
        as.is = TRUE
    )
    GRs <- paste(paste(GRs[, 1], GRs[, 2], sep = ":"), GRs[, 3], sep = "-")
    paste("--intervals", GRs)
}

#' Run GATK HaplotypeCaller
#' @description Internal helper for gatkCalleR(). Executes the GATK
#'     HaplotypeCaller and checks the exit code.
#' @param gatk_path character string with the GATK executable path.
#' @param args character vector with GATK arguments.
#' @param path1 character string with the BAM directory.
#' @param bamFile character string with the BAM file name.
#' @param SE logical, save standard error?
#' @return invisible(NULL). Stops on non-zero exit code.
#' @author Ilia M. Ilyushonak, Vasily V. Grinev.
#' @keywords internal
#' @noRd
.runHaplotypeCaller <- function(gatk_path, args, path1, bamFile, SE) {
    se <- ""
    if (SE == TRUE) {
        se <- paste(
            gsub(pattern = ".bam", replacement = "", x = bamFile),
            "STDERR.HC.TXT",
            sep = "."
        )
    }
    status <- system2(
        command = gatk_path,
        args = args,
        if (se == "") {
            stderr <- se
        } else {
            stderr <- paste(path1, se, sep = "/")
        }
    )
    if (status != 0) {
        stop(
            "GATK HaplotypeCaller failed with exit code ", status,
            ". Check that gatk_path points to a valid GATK executable ",
            "and that the input files (BAM, FASTA, intervals) are correct."
        )
    }
    invisible(NULL)
}


#' Call of single nucleotide variations using Fisher’s exact test
#' @description Internal wrapper around Rsubread::exactSNP() for SNV
#'     calling. Called from callSNVs(criteria = "fisherRSubread").
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character string giving the name of input BAM file.
#' @param fastaDir character string giving the name (or path to and name) of
#'     folder for storing of reference genome.
#' @param faFile charater string giving the name of file containing reference
#'     genome in format *.FASTA or *.FA.
#' @param baseq integer value, minimal "base quality" for each nucleotide in an
#'     alignment. Default value is 20 that means read nucleotides with quality
#'     scores less than 20 will be excluded from analysis.
#' @param mindepth integer value giving the minimal coverage of location
#'     containing the variation. Default value is 1.
#' @param maxdepth integer value giving the maximal coverage of location
#'     containing the variation. Default value is 1e6. This option is useful
#'     for removing PCR artifacts.
#' @param qvalue numeric value giving the q-value cutoff for variations calling
#'     at sequencing depth 50X. Default value is 12. The q-value is calcuated
#'     as -log10(p), where p is a p-value yielded from the Fisher’s exact test.
#'     The function exactSNP() automatically adjusts the q-value cutoff for
#'     each chromosomal location according to its sequencing depth, based on
#'     this cutoff.
#' @param trim integer value giving the number of nucleotides trimmed off from
#'     each end of the read. Default value is 0.
#' @param threads numeric value giving the number of threads/CPUs used. Default
#'     value is 1.
#' @param outputVCF character string giving the name of output VCF file. NULL
#'     by default, in which case the file name is generated automatically as
#'     <bam_basename>.fisherCalleR.SNVs.vcf. If not NULL, the file name is
#'     <bam_basename>.<outputVCF>.vcf.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return file of VCF format is written to the BAM directory. The function
#'     is called for its side effect and returns \code{NULL} invisibly.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
#' @importFrom Rsubread exactSNP
fisherCalleR <- function(bamDir = NULL,
    bamFile,
    fastaDir = NULL,
    faFile,
    baseq = 20,
    mindepth = 1, maxdepth = 1e6,
    qvalue = 12,
    trim = 0,
    threads = 1,
    outputVCF = NULL,
    workDir = NULL) {
    ### Setting the working directory.
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    ### Full path to the BAM file.
    path1 <- .resolveBamPath(bamDir, workDir)
    ### Full path to the reference genome.
    path2 <- .resolveFastaPath(fastaDir, workDir)
    ### Setting the BAM file.
    bam <- paste(path1, bamFile, sep = "/")
    ### Setting the FASTA file.
    fa <- paste(path2, faFile, sep = "/")
    ### Setting the output VCF file.
    if (is.null(x = outputVCF)) {
        vcf <- paste0(sub(".bam$", "", bam), ".fisherCalleR.SNVs.vcf")
    } else {
        vcf <- paste0(sub(".bam$", "", bam), ".", outputVCF, ".vcf")
    }
    ### Detection of SNVs.
    SNVs <- exactSNP(
        readFile = bam,
        isBAM = TRUE,
        refGenomeFile = fa,
        minBaseQuality = baseq,
        minReads = mindepth,
        maxReads = maxdepth,
        qvalueCutoff = qvalue,
        nTrimmedBases = trim,
        nthreads = threads,
        outputFile = vcf
    )
    ### Returning nothing: the VCF file is the only output.
    invisible(NULL)
}


#' Call single nucleotide variations using different statistical tests
#' @description This function is a wrapper that allows to call and control
#'     various low-level functions for identification of single nucleotide
#'     variations from NGS data.
#' @param fr_data either a name of CSV file containing pileup statistics (as
#'     returned by collectBasesPileup() function) or an object of class
#'     data.frame with frequancy table.
#' @param criteria character string specifying which test to be performed.
#'     Only one criterion per call is allowed. The following values are
#'     possible: "binomial", "counting", "entropy" (by default), "fisher",
#'     "fisherRSubread", "gatk", "ml" and "poisson".
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory for storing of generated VCF file(-s). NULL by default, which
#'     means the current working directory will be used instead of the
#'     specified BAM directory.
#' @param outputVCF character string giving the name of output VCF file.
#' @param bamDir character string giving the name (or path to and name) of
#'     directory with BAM file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified BAM directory.
#' @param bamFile character string giving the name of input BAM file.
#' @param fastaDir character string giving the name (or path to and name) of
#'     folder for storing of reference genome (with .fai indices if GATK).
#' @param faFile charater string giving the name of file containing reference
#'     genome in format *.FASTA or *.FA.
#' @param baseq integer value, minimal "base quality" for each nucleotide in an
#'     alignment. Default value is 20 that means read nucleotides with quality
#'     scores less than 20 will be excluded from analysis.
#' @param mindepth integer value giving the minimal coverage of location
#'     containing the variation. Default value is 1.
#' @param maxdepth integer value giving the maximal coverage of location
#'     containing the variation. Default value is 1e6. This option is useful
#'     for removing PCR artifacts.
#' @param qvalue numeric value giving the q-value cutoff for variations calling
#'     at sequencing depth 50X. Default value is 12. The q-value is calcuated
#'     as -log10(p), where p is a p-value yielded from the Fisher’s exact test.
#'     The function exactSNP() automatically adjusts the q-value cutoff for
#'     each chromosomal location according to its sequencing depth, based on
#'     this cutoff.
#' @param trim integer value giving the number of nucleotides trimmed off from
#'     each end of the read. Default value is 0.
#' @param threads numeric value giving the number of threads/CPUs used. Default
#'     value is 1.
#' @param postfix character string that will be added to the output VCF file
#'     name as postfix. NA by default at which the postfix "GATK based SNVs"
#'     will be used. Only with GATK.
#' @param intervals character string specifying the name of the tab-delimited
#'     TXT file containing coordinates of genomic location(-s) of interest. The
#'     default value is NULL. If so, this file must contains the following four
#'     fields:
#'     i) seqnames - the name of chromosome or scaffold with prefix "chr";
#'     ii) start   - start coordinate of genomic interval of interest;
#'     iii) end    - end coordinate of genomic interval of interest;
#'     iv) strand  - (optionally) strand information about genomic interval
#'                   of interest.
#'     TXT file with genomic location(-s) must be in working directory.
#' @param pcr_model GATK argument --pcr-indel-model. Possible values: NONE,
#'     AGRESSIVE, HOSTILE, CONSERVATIVE. See HaplotypeCaller
#'     at https://gatk.broadinstitute.org for further details.
#' @param soft_clip GATK argument --dont-use-soft-clipped. TRUE by default.
#' @param ERC GATK ERC mode. "NONE" by default. Possible values: NONE, GVCF or
#'     BP_RESOLUTION. See HaplotypeCaller at https://gatk.broadinstitute.org
#'     for further details.
#' @param cores number of CPU cores to use for Hidden Markov Models. If "MAX"
#'     (by default), all available cores will be detected and used.
#' @param SE if TRUE (by default), standart error determined by GATK will be
#'     writen in file STDERR.HC.TXT.
#' @param gatk_path character string giving the path to internal or external
#'     (system) GATK executive script. This path can be obtained by running
#'     Find.Gatk() or settled manually. Default value is "gatk", i.e. function
#'     will uses GATK systemly from $PATH.
#' @param workDir character string giving the path to and name of work
#'     directory. NULL by default that means the current working directory.
#' @param model character string specifying which machine learning model to
#'     use. Available options: "log_reg", "rand_forest", "rpart", "SVM",
#'     "xgboost". Default is "rpart". Only with criterion "ml".
#' @param pvalue logical indicating whether to return class probabilities
#'     instead of class labels. Default is FALSE. Only with criterion "ml".
#' @return VCF file(-s) with detected single nucleotide variations are
#'     written to \code{vcfDir}. The function is called for its side effect
#'     and returns \code{NULL} invisibly.
#' @details
#'     The function supports two groups of criteria:
#'     \itemize{
#'       \item Table-based criteria (\code{"binomial"}, \code{"counting"},
#'             \code{"entropy"}, \code{"fisher"}, \code{"poisson"},
#'             \code{"ml"}) require \code{fr_data} --- a frequency table
#'             produced by \code{collectBasesPileup()}.
#'       \item BAM-based criteria (\code{"fisherRSubread"}, \code{"gatk"})
#'             require \code{bamFile} and work directly with the BAM file.
#'     }
#' @author Dzianis D. Sarnatski, Mikalai M. Yatskou, Vasily V. Grinev
#' @examples
#' res <- callSNVs(
#'     fr_data = "test_seq.PileupResults.csv",
#'     criteria = "fisher",
#'     vcfDir = "Files_VCF",
#'     outputVCF = "test_seq.fisher",
#'     workDir = "D:/Vasily Grinev"
#' )
#' @export
#' @importFrom data.table data.table
callSNVs <- function(
        fr_data,
        criteria = "entropy",
        vcfDir = NULL,
        outputVCF,
        bamDir = NULL,
        bamFile = NA,
        fastaDir = NULL,
        faFile = NA,
        baseq = 20,
        mindepth = 1,
        maxdepth = 1e6,
        qvalue = 12,
        trim = 0,
        threads = 1,
        postfix = NA,
        intervals = NULL,
        pcr_model = NA,
        soft_clip = TRUE,
        ERC = "NONE",
        cores = "MAX",
        SE = TRUE,
        gatk_path = NULL,
        workDir = NULL,
        model = "rpart",
        pvalue = FALSE) {
    ### Validation of the input criteria and arguments.
    criteria <- .validateCallerInputs(
        criteria, fr_data, missing(fr_data), bamFile
    )
    ### Setting the working directory.
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    ### Reading the frequency table if needed.
    fr_table <- .readCallerFrequencyTable(criteria, fr_data, workDir)
    ### Dispatch to the low-level function.
    invisible(.dispatchCaller(
        criteria, fr_data, fr_table, vcfDir, outputVCF,
        bamDir, bamFile, fastaDir, faFile,
        baseq, mindepth, maxdepth, qvalue, trim, threads,
        postfix, intervals, pcr_model, soft_clip, ERC, cores, SE,
        gatk_path, workDir, model, pvalue
    ))
    ### Returning nothing.
    invisible(NULL)
}

#' Validate caller criteria and required arguments
#' @description Internal helper for callSNVs(). Validates the criteria
#'     against SUPPORTED_CALLERS and checks that the required arguments
#'     are supplied for the chosen caller group.
#' @param criteria character string with the chosen criterion.
#' @param fr_data original fr_data argument.
#' @param fr_data_missing logical, whether fr_data is missing.
#' @param bamFile character string with the BAM file name, or NA.
#' @return character string with the validated criterion.
#' @author Dzianis D. Sarnatski, Mikalai M. Yatskou, Vasily V. Grinev.
#' @keywords internal
#' @noRd
.validateCallerInputs <- function(criteria, fr_data, fr_data_missing,
    bamFile) {
    criteria <- match.arg(arg = criteria, choices = SUPPORTED_CALLERS)
    if (criteria %in% SUPPORTED_TABLE_CALLERS && fr_data_missing) {
        stop(
            "fr_data is required for criteria = '", criteria, "'. ",
            "Use collectBasesPileup() to generate it from a BAM file."
        )
    }
    if (criteria %in% SUPPORTED_BAM_CALLERS && is.na(bamFile)) {
        stop("bamFile is required for criteria = '", criteria, "'.")
    }
    criteria
}

#' Read the frequency table for table-based callers
#' @description Internal helper for callSNVs(). Reads the frequency table
#'     when a table-based criterion is chosen; otherwise returns NULL.
#' @param criteria character string with the chosen criterion.
#' @param fr_data original fr_data argument.
#' @param workDir character string with the working directory.
#' @return data.frame with the frequency table, or NULL.
#' @author Dzianis D. Sarnatski, Mikalai M. Yatskou, Vasily V. Grinev.
#' @keywords internal
#' @noRd
.readCallerFrequencyTable <- function(criteria, fr_data, workDir) {
    if (criteria %in% SUPPORTED_TABLE_CALLERS) {
        return(.readFrequencyTable(fr_data, workDir))
    }
    NULL
}

#' Read frequency table from CSV or return data.frame
#' @description Internal helper for callSNVs(). Reads the frequency table
#'     from a CSV file or returns the data.frame as-is.
#' @param fr_data character string with the file name or data.frame.
#' @param workDir character string with the working directory.
#' @return data.frame with the frequency table.
#' @author Dzianis D. Sarnatski, Mikalai M. Yatskou, Vasily V. Grinev.
#' @keywords internal
#' @noRd
.readFrequencyTable <- function(fr_data, workDir) {
    if (is.object(x = fr_data)) {
        return(fr_data)
    }
    path <- paste(workDir, fr_data, sep = "/")
    frt <- tools::file_ext(x = path)
    if (!frt %in% "csv") {
        stop("Invalid format of CSV file")
    }
    if (!file.exists(path)) {
        stop("File not found:", path)
    }
    read.csv(file = path, header = TRUE)
}

#' Dispatch SNV calling to the appropriate low-level function
#' @description Internal helper for callSNVs(). Dispatches the SNV calling
#'     to the low-level function corresponding to the chosen criterion.
#' @param criteria character string with the chosen criterion.
#' @param fr_data original fr_data argument.
#' @param fr_table frequency table, or NULL.
#' @param vcfDir character string with the VCF directory.
#' @param outputVCF character string with the output VCF name.
#' @param bamDir character string with the BAM directory.
#' @param bamFile character string with the BAM file name.
#' @param fastaDir character string with the FASTA directory.
#' @param faFile character string with the FASTA file name.
#' @param baseq integer, minimal base quality.
#' @param mindepth integer, minimal coverage.
#' @param maxdepth integer, maximal coverage.
#' @param qvalue numeric, q-value cutoff.
#' @param trim integer, number of trimmed bases.
#' @param threads integer, number of threads.
#' @param postfix character string with the output postfix.
#' @param intervals character string with the intervals file, or NULL.
#' @param pcr_model character string with the PCR indel model.
#' @param soft_clip logical, use soft-clipped bases?
#' @param ERC character string with the ERC mode.
#' @param cores number or "MAX".
#' @param SE logical, save standard error?
#' @param gatk_path character string with the GATK executable.
#' @param workDir character string with the working directory.
#' @param model character string with the ML model name.
#' @param pvalue logical, return probabilities?
#' @return result of the low-level function.
#' @author Dzianis D. Sarnatski, Mikalai M. Yatskou, Vasily V. Grinev.
#' @keywords internal
#' @noRd
.dispatchCaller <- function(criteria, fr_data, fr_table, vcfDir, outputVCF,
    bamDir, bamFile, fastaDir, faFile,
    baseq, mindepth, maxdepth, qvalue, trim, threads,
    postfix, intervals, pcr_model, soft_clip, ERC,
    cores, SE, gatk_path, workDir, model, pvalue) {
    switch(criteria,
        binomial = ,
        counting = ,
        entropy = ,
        fisher = ,
        poisson = .dispatchTableCaller(
            criteria, fr_table, vcfDir, outputVCF, workDir
        ),
        fisherRSubread = .dispatchFisherRSubread(
            bamDir, bamFile, fastaDir, faFile,
            baseq, mindepth, maxdepth, qvalue, trim, threads,
            outputVCF, workDir
        ),
        gatk = .dispatchGatk(
            bamDir, bamFile, postfix, fastaDir, faFile, intervals,
            pcr_model, soft_clip, ERC, cores, SE, gatk_path, workDir
        ),
        ml = .dispatchMl(
            fr_table, model, pvalue, vcfDir, outputVCF, workDir
        )
    )
}

#' Dispatch table-based callers
#' @description Internal helper for .dispatchCaller(). Runs the chosen
#'     statistical test, appends its result to the frequency table and
#'     generates the VCF file.
#' @param criteria character string with the test name.
#' @param fr_table data.frame with the frequency table.
#' @param vcfDir character string with the VCF directory.
#' @param outputVCF character string with the output VCF name.
#' @param workDir character string with the working directory.
#' @return invisible(NULL).
#' @author Dzianis D. Sarnatski, Mikalai M. Yatskou, Vasily V. Grinev.
#' @keywords internal
#' @noRd
.dispatchTableCaller <- function(criteria, fr_table, vcfDir, outputVCF,
    workDir) {
    result <- .runStatisticalTest(criteria, fr_table)
    input <- data.table(cbind(fr_table, result))
    colnames(input)[ncol(input)] <- criteria
    generateVCF(
        input = input,
        version = "4.0",
        vcfDir = vcfDir,
        outputVCF = outputVCF,
        workDir = workDir
    )
    invisible(NULL)
}

#' Dispatch fisherRSubread caller
#' @description Internal helper for .dispatchCaller(). Wraps the
#'     fisherCalleR() call with its parameters.
#' @param bamDir character string with the BAM directory.
#' @param bamFile character string with the BAM file name.
#' @param fastaDir character string with the FASTA directory.
#' @param faFile character string with the FASTA file name.
#' @param baseq integer, minimal base quality.
#' @param mindepth integer, minimal coverage.
#' @param maxdepth integer, maximal coverage.
#' @param qvalue numeric, q-value cutoff.
#' @param trim integer, number of trimmed bases.
#' @param threads integer, number of threads.
#' @param outputVCF character string with the output VCF name.
#' @param workDir character string with the working directory.
#' @return invisible(NULL).
#' @author Dzianis D. Sarnatski, Mikalai M. Yatskou, Vasily V. Grinev.
#' @keywords internal
#' @noRd
.dispatchFisherRSubread <- function(bamDir, bamFile, fastaDir, faFile,
    baseq, mindepth, maxdepth, qvalue,
    trim, threads, outputVCF, workDir) {
    fisherCalleR(
        bamDir = bamDir,
        bamFile = bamFile,
        fastaDir = fastaDir,
        faFile = faFile,
        baseq = baseq,
        mindepth = mindepth,
        maxdepth = maxdepth,
        qvalue = qvalue,
        trim = trim,
        threads = threads,
        outputVCF = outputVCF,
        workDir = workDir
    )
    invisible(NULL)
}

#' Dispatch GATK caller
#' @description Internal helper for .dispatchCaller(). Wraps the
#'     gatkCalleR() call with its parameters.
#' @param bamDir character string with the BAM directory.
#' @param bamFile character string with the BAM file name.
#' @param postfix character string with the output postfix.
#' @param fastaDir character string with the FASTA directory.
#' @param faFile character string with the FASTA file name.
#' @param intervals character string with the intervals file, or NULL.
#' @param pcr_model character string with the PCR indel model.
#' @param soft_clip logical, use soft-clipped bases?
#' @param ERC character string with the ERC mode.
#' @param cores number or "MAX".
#' @param SE logical, save standard error?
#' @param gatk_path character string with the GATK executable.
#' @param workDir character string with the working directory.
#' @return invisible(NULL).
#' @author Dzianis D. Sarnatski, Mikalai M. Yatskou, Vasily V. Grinev.
#' @keywords internal
#' @noRd
.dispatchGatk <- function(bamDir, bamFile, postfix, fastaDir, faFile,
    intervals, pcr_model, soft_clip, ERC, cores, SE,
    gatk_path, workDir) {
    gatkCalleR(
        bamDir = bamDir,
        bamFile = bamFile,
        postfix = postfix,
        fastaDir = fastaDir,
        faFile = faFile,
        intervals = intervals,
        pcr_model = pcr_model,
        soft_clip = soft_clip,
        ERC = ERC,
        cores = cores,
        SE = SE,
        gatk_path = gatk_path,
        workDir = workDir
    )
    invisible(NULL)
}

#' Dispatch ML caller
#' @description Internal helper for .dispatchCaller(). Wraps the
#'     mlCalleR() call with its parameters.
#' @param fr_table data.frame with the frequency table.
#' @param model character string with the ML model name.
#' @param pvalue logical, return probabilities?
#' @param vcfDir character string with the VCF directory.
#' @param outputVCF character string with the output VCF name.
#' @param workDir character string with the working directory.
#' @return invisible(NULL).
#' @author Dzianis D. Sarnatski, Mikalai M. Yatskou, Vasily V. Grinev.
#' @keywords internal
#' @noRd
.dispatchMl <- function(fr_table, model, pvalue, vcfDir, outputVCF,
    workDir) {
    mlCalleR(
        fr_table = fr_table,
        model = model,
        pvalue = pvalue,
        vcfDir = vcfDir,
        outputVCF = outputVCF,
        workDir = workDir
    )
    invisible(NULL)
}

#' Run a statistical test on the frequency table
#' @description Internal helper for callSNVs(). Dispatches to the appropriate
#'     C++ statistical test.
#' @param criteria character string with the test name.
#' @param fr_table data.frame with the frequency table.
#' @return numeric vector with test results.
#' @author Dzianis D. Sarnatski, Mikalai M. Yatskou, Vasily V. Grinev.
#' @keywords internal
#' @noRd
.runStatisticalTest <- function(criteria, fr_table) {
    if (criteria == "entropy") {
        return(snp_entropy_cpp_parallel(fr_table, p_value = TRUE))
    }
    tests <- list(
        counting = snp_counting_cpp_parallel,
        binomial = snp_binomial_cpp_parallel,
        fisher = snp_fisher_cpp_parallel,
        poisson = snp_poisson_cpp_parallel
    )
    tests[[criteria]](fr_table)
}


#' Call single nucleotide variations using machine learning models
#' @description Internal wrapper around pre-trained machine learning models
#'     for SNV calling. Called from callSNVs(criteria = "ml").
#' @param fr_table data.frame with frequency table.
#' @param model character string specifying which machine learning model to use.
#'     Available options: "log_reg", "rand_forest", "rpart", "SVM", "xgboost".
#'     Default is "rpart".
#' @param pvalue logical indicating whether to return class probabilities
#'     instead of class labels (default: FALSE). If TRUE, returns probability
#'     of belonging to positive class (1).
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory for storing of generated VCF file(-s).
#' @param outputVCF character string giving the name of output VCF file.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return VCF file is written to \code{vcfDir}. The function is called for
#'     its side effect and returns \code{NULL} invisibly.
#' @author Dzianis D. Sarnatski
#' @keywords internal
#' @noRd
#' @importFrom data.table data.table
#' @import randomForest
#' @import LiblineaR
#' @import rpart
mlCalleR <- function(fr_table,
    model = "rpart",
    pvalue = FALSE,
    vcfDir = NULL,
    outputVCF,
    workDir = NULL) {
    ### Validation of the input model.
    if (!model %in% SUPPORTED_ML_MODELS) {
        stop(paste(
            "Invalid model specified:", model,
            ". Available options:",
            paste(SUPPORTED_ML_MODELS, collapse = ", ")
        ))
    }
    ### Setting the working directory.
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    ### Validate input data structure.
    required_cols <- c("position", "reference", "A", "C", "G", "T")
    if (!all(required_cols %in% colnames(fr_table))) {
        stop(paste(
            "Input data must contain columns:",
            paste(required_cols, collapse = ", ")
        ))
    }
    ### Prepare features.
    prepared <- .prepareMLFeatures(fr_table)
    ### Load model.
    model.obj <- .loadMLModel(model)
    ### Make predictions.
    prediction <- .predictML(model.obj, model, prepared$data.features, pvalue)
    ### Add predictions to the original data.
    results <- data.table(cbind(
        prepared$fr_table_filt,
        ML_Prediction = prediction
    ))
    ### Generate VCF file.
    generateVCF(
        input = results,
        version = "4.0",
        vcfDir = vcfDir,
        outputVCF = outputVCF,
        workDir = workDir
    )
    ### Returning nothing.
    invisible(NULL)
}

#' Prepare features for ML prediction
#' @description Internal helper for mlCalleR(). Extracts the features from the
#'     frequency table, normalizes them, and filters out non-finite rows.
#' @param fr_table data.frame with the frequency table.
#' @return list with data.features and fr_table_filt.
#' @author Dzianis D. Sarnatski.
#' @keywords internal
#' @noRd
.prepareMLFeatures <- function(fr_table) {
    features <- c("x1", "x2", "x3", "x4")
    data.features <- vectorise_data_cpp_parallel(fr_table, features = features)
    total <- rowSums(data.features)
    target_cols <- c("x2", "x3", "x4")
    existing_cols <- target_cols[target_cols %in% names(data.features)]
    if (length(existing_cols) > 0) {
        x_sum <- rowSums(
            data.features[, existing_cols, drop = FALSE],
            na.rm = TRUE
        )
    } else {
        stop("Non-reference columns (x2, x3, x4) are missing in data.features.")
    }
    data.features <- data.frame(
        x_sum_norm = x_sum / total,
        total = total
    )
    finite_rows <- is.finite(data.features$x_sum_norm) &
        is.finite(data.features$total)
    list(
        data.features = data.features[finite_rows, , drop = FALSE],
        fr_table_filt = fr_table[finite_rows, , drop = FALSE]
    )
}

#' Load pre-trained ML model
#' @description Internal helper for mlCalleR(). Loads the pre-trained model
#'     from inst/extdata/.
#' @param model character string with the model name.
#' @return loaded model object.
#' @author Dzianis D. Sarnatski.
#' @keywords internal
#' @noRd
.loadMLModel <- function(model) {
    if (model == "xgboost") {
        model.path <- system.file("extdata", paste0(model, ".json"),
            package = "huGSVcalleR"
        )
        if (!nzchar(model.path) || !file.exists(model.path)) {
            stop(
                "Model file not found: ", model.path,
                ". Expected a .json file in inst/extdata/."
            )
        }
        xgboost::xgb.load(model.path)
    } else {
        model.path <- system.file("extdata", paste0(model, ".rds"),
            package = "huGSVcalleR"
        )
        if (!nzchar(model.path) || !file.exists(model.path)) {
            stop(
                "Model file not found: ", model.path,
                ". Expected a .rds file in inst/extdata/."
            )
        }
        readRDS(model.path)
    }
}

#' Make ML predictions
#' @description Internal helper for mlCalleR(). Dispatches to the appropriate
#'     prediction function based on the model type.
#' @param model.obj loaded model object.
#' @param model character string with the model name.
#' @param data.features data.frame with the features.
#' @param pvalue logical, return probabilities?
#' @return prediction vector.
#' @author Dzianis D. Sarnatski.
#' @keywords internal
#' @noRd
.predictML <- function(model.obj, model, data.features, pvalue) {
    if (model == "xgboost") {
        .predictXgboost(model.obj, data.features, pvalue)
    } else if (model == "rpart") {
        .predictRpart(model.obj, data.features, pvalue)
    } else if (model == "log_reg") {
        .predictLogReg(model.obj, data.features, pvalue)
    } else if (model == "rand_forest") {
        .predictRandForest(model.obj, data.features, pvalue)
    } else if (model == "SVM") {
        .predictSVM(model.obj, data.features, pvalue)
    }
}

#' Predict with xgboost
#' @description Internal helper for mlCalleR(). Makes predictions with an
#'     xgboost model.
#' @param model.obj loaded xgboost model.
#' @param data.features data.frame with the features.
#' @param pvalue logical, return probabilities?
#' @return prediction vector.
#' @author Dzianis D. Sarnatski.
#' @keywords internal
#' @noRd
.predictXgboost <- function(model.obj, data.features, pvalue) {
    pred <- predict(
        model.obj,
        xgboost::xgb.DMatrix(data = as.matrix(data.features))
    )
    if (pvalue) {
        pred
    } else {
        pred |>
            round() |>
            as.factor()
    }
}

#' Predict with rpart
#' @description Internal helper for mlCalleR(). Makes predictions with an
#'     rpart model.
#' @param model.obj loaded rpart model.
#' @param data.features data.frame with the features.
#' @param pvalue logical, return probabilities?
#' @return prediction vector.
#' @author Dzianis D. Sarnatski.
#' @keywords internal
#' @noRd
.predictRpart <- function(model.obj, data.features, pvalue) {
    if (pvalue) {
        pred <- predict(model.obj, newdata = data.features, type = "prob")
        if (is.list(pred)) {
            vapply(pred, function(x) x[2], numeric(1))
        } else if (is.matrix(pred)) {
            pred[, 2]
        } else {
            pred
        }
    } else {
        predict(model.obj, newdata = data.features, type = "class")
    }
}

#' Predict with logistic regression
#' @description Internal helper for mlCalleR(). Makes predictions with a
#'     logistic regression model.
#' @param model.obj loaded glm model.
#' @param data.features data.frame with the features.
#' @param pvalue logical, return probabilities?
#' @return prediction vector.
#' @author Dzianis D. Sarnatski.
#' @keywords internal
#' @noRd
.predictLogReg <- function(model.obj, data.features, pvalue) {
    pred <- predict(
        model.obj,
        newdata = data.features,
        type = "response"
    ) |> as.vector()
    if (pvalue) {
        pred
    } else {
        ifelse(pred > 0.5, 1, 0) |> as.factor()
    }
}

#' Predict with random forest
#' @description Internal helper for mlCalleR(). Makes predictions with a
#'     random forest model.
#' @param model.obj loaded randomForest model.
#' @param data.features data.frame with the features.
#' @param pvalue logical, return probabilities?
#' @return prediction vector.
#' @author Dzianis D. Sarnatski.
#' @keywords internal
#' @noRd
.predictRandForest <- function(model.obj, data.features, pvalue) {
    if (pvalue) {
        pred <- predict(model.obj, newdata = data.features, type = "prob")
        if (is.list(pred)) {
            vapply(pred, function(x) x[2], numeric(1))
        } else if (is.matrix(pred)) {
            pred[, 2]
        } else {
            as.vector(pred)
        }
    } else {
        predict(model.obj, newdata = data.features, type = "class")
    }
}

#' Predict with SVM
#' @description Internal helper for mlCalleR(). Makes predictions with an
#'     SVM model.
#' @param model.obj loaded LiblineaR model.
#' @param data.features data.frame with the features.
#' @param pvalue logical, return probabilities?
#' @return prediction vector.
#' @author Dzianis D. Sarnatski.
#' @keywords internal
#' @noRd
.predictSVM <- function(model.obj, data.features, pvalue) {
    if (pvalue) {
        pred_obj <- predict(
            model.obj,
            newx = as.matrix(data.features),
            proba = TRUE
        )
        if (is.null(pred_obj$probabilities) ||
            ncol(pred_obj$probabilities) < 2) {
            stop(
                "SVM model does not support probability output. ",
                "It must be trained with LiblineaR type 0, 6 or 7 ",
                "(logistic regression) to compute probabilities. ",
                "Use pvalue = FALSE for class labels."
            )
        }
        pred_obj$probabilities[, "1"]
    } else {
        predict(
            model.obj,
            newx = as.matrix(data.features),
            proba = FALSE
        )$predictions
    }
}


#' Generate VCF file based on call table
#' @description Internal helper that converts a call table to VCF. Called
#'     from callSNVs() and mlCalleR().
#' @param input character string giving the name of tab-delimited TXT file or
#'     R object of class data.frame containing call data. File or object must
#'     contains the following fields:
#'     i) position   - location of single nucleotide variant in format
#'                     chrZ:coordinate;
#'     ii) reference - reference variant of nucleotide;
#'     iii) A        - number of redas supporting A in a given position;
#'     iv) C         - number of redas supporting C in a given position;
#'     v) G          - number of redas supporting G in a given position;
#'     vi) T         - number of redas supporting T in a given position;
#'     vii-xi) ...   - one to five columns with statistical test results.
#'     If file, it must be located in working directory.
#' @param version character string specifying the version of VCF format.
#'     Default value is "4.0". It is planned that other versions of the format
#'     will be added in the future.
#' @param vcfDir character string giving the name (or path to and name) of
#'     directory with VCF file(-s). NULL by default, which means the current
#'     working directory will be used instead of the specified VCF directory.
#' @param outputVCF character string giving the name of output VCF file.
#' @param workDir character string giving the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return file(-s) of VCF format is/are written to the current working
#'     directory or specified VCF directory. The function is called for its
#'     side effect and returns \code{NULL} invisibly.
#' @details
#'     For tests that return a p-value (\code{"binomial"}, \code{"poisson"},
#'     \code{"fisher"}, \code{"entropy"}), the QUAL field is written as a
#'     Phred-scaled score \code{-10 * log10(p)}, which is the standard VCF
#'     convention. Values with \code{p = 0} (underflow) are written as 1000.
#'     For \code{"counting"} the QUAL field contains the raw ratio of
#'     reference reads (0--1), because this test does not return a p-value.
#'     For \code{"ml"} the QUAL field contains the probability of the
#'     positive class (0--1).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
#' @importFrom data.table fread
generateVCF <- function(
        input,
        version = "4.0",
        vcfDir = NULL,
        outputVCF,
        workDir = NULL) {
    ### Setting the working directory.
    if (is.null(x = workDir)) {
        workDir <- getwd()
    }
    ### Setting the VCF directory.
    if (is.null(x = vcfDir)) {
        vcfDir <- workDir
    } else {
        vcfDir <- paste(workDir, vcfDir, sep = "/")
    }
    ### Retrieving of the input data.
    if (is.object(x = input)) {
        vcf_pre <- as.data.table(input)
    } else {
        vcf_pre <- fread(file = paste(workDir, input, sep = "/"))
    }
    ### Test columns (everything after the standard 6 columns).
    tests <- colnames(x = vcf_pre)[-seq_len(6)]
    ### Transformation of data into a frame of the appropriate structure.
    VCF <- .prepareVcfTable(vcf_pre, tests)
    if (version == "4.0") {
        ### Writing of VCF file(-s).
        .writeVcfFiles(VCF, tests, vcfDir, outputVCF)
    }
    ### Returning nothing.
    invisible(NULL)
}

#' Prepare VCF table from call data
#' @description Internal helper for generateVCF(). Transforms the call table
#'     into the VCF-like structure with CHROM, POS, ID, REF, ALT, QUAL,
#'     FILTER, INFO columns and per-test result columns.
#' @param vcf_pre data.table with the call data.
#' @return data.frame with the VCF structure.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.prepareVcfTable <- function(vcf_pre, tests) {
    vcf_pre$DP <- rowSums(x = vcf_pre[, c(3, 4, 5, 6)])
    A <- vcf_pre[vcf_pre$A > 0 & vcf_pre$reference != "A", c(-4, -5, -6)]
    A$ALT <- "A"
    A$MM <- A$A
    A <- A[, -3]
    C <- vcf_pre[vcf_pre$C > 0 & vcf_pre$reference != "C", c(-3, -5, -6)]
    C$ALT <- "C"
    C$MM <- C$C
    C <- C[, -3]
    G <- vcf_pre[vcf_pre$G > 0 & vcf_pre$reference != "G", c(-3:-4, -6)]
    G$ALT <- "G"
    G$MM <- G$G
    G <- G[, -3]
    T_sub <- vcf_pre[vcf_pre$T > 0 & vcf_pre$reference != "T", c(-3:-5)]
    T_sub$ALT <- "T"
    T_sub$MM <- T_sub$T
    T_sub <- T_sub[, -3]
    vcf_pre <- rbind(A, C, G, T_sub)
    VCF <- .buildPositionSummary(vcf_pre)
    VCF$CHROM <- gsub(pattern = ":.+", replacement = "", x = VCF$position)
    VCF$POS <- as.numeric(x = gsub(
        pattern = ".+:",
        replacement = "",
        x = VCF$position
    ))
    VCF$ID <- paste(
        paste(VCF$position, VCF$reference, sep = "_"),
        VCF$ALT,
        sep = "/"
    )
    VCF$QUAL <- "."
    VCF$FILTER <- "."
    VCF$INFO <- paste(
        paste0("DP=", VCF$DP),
        paste0("MMsum=", VCF$MMsum),
        paste0("MM=", VCF$MM),
        sep = ";"
    )
    VCF$REF <- VCF$reference
    VCF <- data.frame(VCF)
    ### Keep the test columns along with the standard VCF columns.
    VCF[, c(
        "CHROM", "POS", "ID", "REF", "ALT",
        "QUAL", "FILTER", "INFO", tests
    )]
}

#' Build per-position ALT/MMsum/MM summary
#' @description Internal helper for .prepareVcfTable(). Collapses the long
#'     format table into one row per position, aggregating ALT alleles,
#'     total supporting reads (MMsum) and per-allele supporting reads (MM).
#' @param vcf_pre data.table in long format (one row per ALT allele).
#' @return data.frame with position, reference, DP, ALT, MMsum and MM.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.buildPositionSummary <- function(vcf_pre) {
    VCF <- unique(as.data.frame(vcf_pre)[, seq_len(ncol(x = vcf_pre) - 2)])
    VCF <- VCF[order(x = VCF$position), ]
    ALT <- vcf_pre[, paste(ALT, collapse = ","), by = position]
    ALT <- ALT[order(x = ALT$position), ]
    VCF$ALT <- ALT$V1
    MMsum <- vcf_pre[, sum(MM), by = position]
    MMsum <- MMsum[order(x = MMsum$position), ]
    VCF$MMsum <- MMsum$V1
    MM <- vcf_pre[, paste(MM, collapse = ","), by = position]
    MM <- MM[order(x = MM$position), ]
    VCF$MM <- MM$V1
    VCF
}

#' Write VCF files for each test
#' @description Internal helper for generateVCF(). Writes one VCF file per
#'     test, applying the appropriate QUAL transformation.
#' @param VCF data.frame with the VCF structure.
#' @param tests character vector with test names.
#' @param vcfDir character string with the VCF directory.
#' @param outputVCF character string with the output VCF name.
#' @return invisible(NULL).
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.writeVcfFiles <- function(VCF, tests, vcfDir, outputVCF) {
    meta <- .prepareVcfMeta()
    pvalue_tests <- c("binomial", "poisson", "fisher", "entropy")
    for (i in seq_along(tests)) {
        VCF$QUAL <- .computeQualValues(VCF, tests[i], pvalue_tests)
        vcf_name <- paste(
            paste(vcfDir,
                paste(outputVCF, tests[i], sep = "."),
                sep = "/"
            ),
            "vcf",
            sep = "."
        )
        meta_file <- file(description = vcf_name, open = "wt")
        writeLines(text = meta, con = meta_file)
        close(con = meta_file)
        write.table(
            x = VCF[, seq_len(8)],
            file = vcf_name,
            sep = "\t",
            quote = FALSE,
            col.names = FALSE,
            row.names = FALSE,
            append = TRUE
        )
    }
    invisible(NULL)
}

#' Compute QUAL values for a given test
#' @description Internal helper for generateVCF(). Applies the appropriate
#'     transformation to the QUAL column: Phred scaling for p-value tests,
#'     raw values for counting and ML.
#' @param VCF data.frame with the VCF structure.
#' @param test character string with the test name.
#' @param pvalue_tests character vector with the names of p-value tests.
#' @return numeric or character vector with QUAL values.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.computeQualValues <- function(VCF, test, pvalue_tests) {
    qual_values <- VCF[[test]]
    if (is.numeric(qual_values) && test %in% pvalue_tests) {
        qual_values <- pmin(pmax(qual_values, 0), 1)
        qual_values <- -10 * log10(qual_values)
        qual_values[is.infinite(qual_values)] <- 1000
        round(x = qual_values, digits = 4)
    } else if (is.numeric(qual_values)) {
        round(x = qual_values, digits = 4)
    } else {
        as.character(qual_values)
    }
}

#' Prepare VCF meta-information
#' @description Internal helper for generateVCF(). Builds the VCF meta-header
#'     with the INFO field definitions.
#' @return character vector with VCF meta-lines.
#' @author Vasily V. Grinev.
#' @keywords internal
#' @noRd
.prepareVcfMeta <- function() {
    c(
        "##fileformat=VCFv4.0",
        paste0(
            "##INFO=<ID=DP,Number=1,Type=Integer,",
            "Description=\"Read Depth\">"
        ),
        paste0(
            "##INFO=<ID=MMsum,Number=1,Type=Integer,",
            "Description=\"Number of supporting reads",
            " for alternative allele(-s)\">"
        ),
        paste0(
            "##INFO=<ID=MM,Number=1,Type=String,",
            "Description=\"Number of supporting reads",
            " for each alternative allele\">"
        ),
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
    )
}
