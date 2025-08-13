#' Collect and present alignment quality assessment results
#' @description A high-level function for consolidating and presenting the
#'     results of evaluation of the alignment quality metrics for aligned
#'     DAN-Seq and/or RNA-Seq short reads.
#' @param x a character string specifying the name of object developed by the
#'     function assessQAlignedReads().
#' @param output a character string that specifying a postfix for naming all
#'     output report files. The default value is NULL, at which the name of
#'     sample will be used for naming all output files.
#' @param workDir character string specifying the path to and name of working
#'     directory. NULL by default that means the current working directory.
#' @return a series of text and graphic files.
#' @author Vasily V. Grinev
#' @examples
#' report <- reportQAlResults(x=qal,
#'                            output=NULL,
#'                            workDir="D:/Vasily Grinev")
#' @export
#' Last updated: August 13, 2025.

reportQAlResults <- function(x,
                             output=NULL,
                             workDir=NULL){
    ### Loading the required packages.
    #   This code was successfully tested with packages grid, gridExtra v.2.3
    #   and vioplot v.0.5.0.
    suppressMessages(expr=library(package=grid))
    suppressMessages(expr=library(package=gridExtra))
    suppressMessages(expr=library(package=vioplot))
    ### Textual report.
    ##  Template of textual table.
    tab <- data.frame(Parameter=c("Sample name",
                                  "File name",
                                  "Number of reads, total",
                                  "Number of reads, unique",
                                  "Number of reads, duplicated",
                                  "Number of reads, mapped",
                                  "Number of reads, unmapped",
                                  "Number of reads, paired",
                                  "Number of reads, properly paired",
                                  "Number of reads, not properly paired",
                                  "Number of reads, plus strand",
                                  "Number of reads, minus strand",
                                  "Number of reads, singletons",
                                  "Secondary alignments",
                                  "Supplementary (chimeric) alignments",
                                  "Mapping quality score, min",
                                  "Mapping quality score, 1st quartile",
                                  "Mapping quality score, median",
                                  "Mapping quality score, 3rd quartile",
                                  "Mapping quality score, max",
                                  "Template length, min",
                                  "Template length, 1st quartile",
                                  "Template length, median",
                                  "Template length, 3rd quartile",
                                  "Template length, max"),
                      Value=0) 
    ##  Values of textual table.
    tab[1, 2] <- x$SN
    tab[2, 2] <- x$FN
    tab[3, 2] <- x$TNRs
    tab[4, 2] <- x$TNRs - x$FLAGs["isDuplicate"][[1]]
    tab[5, 2] <- x$FLAGs["isDuplicate"][[1]]
    tab[6, 2] <- x$TNRs - x$FLAGs["isUnmappedQuery"][[1]]
    tab[7, 2] <- x$FLAGs["isUnmappedQuery"][[1]]
    tab[8, 2] <- x$FLAGs["isPaired"][[1]]
    tab[9, 2] <- x$FLAGs["isProperPair"][[1]]
    tab[10, 2] <- x$FLAGs["isPaired"][[1]] - x$FLAGs["isProperPair"][[1]]
    tab[11, 2] <- x$TNRs - x$FLAGs["isUnmappedQuery"][[1]] -
                  x$FLAGs["isMinusStrand"][[1]]
    tab[12, 2] <- x$FLAGs["isMinusStrand"][[1]]
    tab[13, 2] <- x$FLAGs["hasUnmappedMate"][[1]]
    tab[14, 2] <- x$FLAGs["isSecondaryAlignment"][[1]]
    tab[15, 2] <- x$FLAGs["isSupplementaryAlignment"][[1]]
    tab[16:20, 2] <- quantile(x=x$MAPQs, probs=c(0.1, 0.25, 0.5, 0.75, 0.9))
    tab[21:25, 2] <- quantile(x=x$TLENs, probs=c(0.1, 0.25, 0.5, 0.75, 0.9))
    ##  Saving of textual table as *.PDF file.
    if (is.null(x=output)){
        file_output <- paste(sprintf("Mapping quality metrics of %s", x$SN),
                             "pdf",
                             sep=".")
    }else{
        file_output <- paste(sprintf("Mapping quality metrics of %s", output),
                             "pdf",
                             sep=".")
    }
    pdf(file=paste(workDir, file_output, sep="/"), height=nrow(x=tab)/3)
    theme=ttheme_minimal(core=list(fg_params=list(hjust=0,
                                                  x=0.03,
                                                  fontsize=10),
                                   bg_params=list(fill=rep(c("azure1",
                                                             "azure3"),
                                                           each = 1),
                                                           col="black")),
                         colhead=list(fg_params=list(col="white"),
                                      bg_params=list(fill="darkcyan",
                                                     col="black")))
    table <- tableGrob(d=tab, rows=NULL, theme=theme)
    title1 <- textGrob(label="Mapping quality. Summary statistics",
                       just="centre",
                       vjust=-32.4,
                       gp=gpar(col="black",
                               font=2,
                               fontsize=12))
    title2 <- textGrob(label=sprintf("Sample: %s", x$SN),
                       just="centre",
                       vjust=-33.5,
                       gp=gpar(col="black",
                               fontsize=11))
    title3 <- textGrob(label=sprintf(paste("Date & Time:",
                                           format(x=Sys.time(), "%Y-%m-%d %X"),
                                           sep=" ")),
                       just="centre",
                       vjust=-34.9,
                       gp=gpar(col="black",
                               fontsize=10))
    note <- "Generated by assessQAlignedReads() & reportQAlResults()"
    note <- textGrob(label=sprintf(note),
                     just="centre",
                     vjust=47,
                     gp=gpar(col="black",
                             fontsize=8))
    gt <- gTree(children=gList(table, title1, title2, title3, note))
    grid.draw(gt)
    suppressMessages(expr=dev.off())
    ### Graphical reports.
    ##  Plots of mapping qualities.
    if (is.null(x=output)){
        file_output <- paste(sprintf("Mapping quality plots of %s", x$SN),
                             "pdf",
                             sep=".")
    }else{
        file_output <- paste(sprintf("Mapping quality plots of %s", output),
                             "pdf",
                             sep=".")
    }
    pdf(file=paste(workDir, file_output, sep="/"),
        onefile=TRUE,
        paper="a4",
        width=8.27, height=11.69)
    par(mfrow=c(3, 2),
        omi=c(0.79, 0.79, 1, 0.79))
    hist(x=x$MAPQs, axes=FALSE, ann=FALSE, border=FALSE,
         xlim=c(0, 40),
         col=rgb(red=0, green=158, blue=115, maxColorValue=255))
    axis(side=1,
         at=c(0, 10, 20, 30, 40),
         labels=c(0, 10, 20, 30, 40),
         las=1, lwd=0.8)
    axis(side=2,
         at=c(0.0, max(x=table(x$MAPQs))/3 * 1:3),
         labels=formatC(x=c(0.0, max(x=table(x$MAPQs))/3 * 1:3),
                        format="e", digits=2),
         lwd=0.8)
    title(main="Distribution of mapping qualities", cex.main=1.5)
    title(xlab="Mapping quality scores", line=2.5, cex.lab=1.2)
    title(ylab="Number of reads", line=3, cex.lab=1.2)
    per <- data.frame(threshold=0:40, percent=0)
    per[1, 2] <- 100
    L <- length(x=x$MAPQs)
    for (i in 1:40){
        per[i + 1, 2] <- length(x=x$MAPQs[x$MAPQs >= i])/L * 100
    }
    plot(per, type="l", bty="n",
         xlim=c(0, 40), ylim=c(0, 100),
         lwd=1,
         col=rgb(red=0, green=114, blue=178, maxColorValue=255),
         ann=FALSE, axes=FALSE)
    axis(side=1,
         at=c(0, 10, 20, 30, 40),
         labels=c(0, 10, 20, 30, 40),
         las=1, lwd=0.8)
    axis(side=2,
         at=c(0, 20, 40, 60, 80, 100),
         labels=c(0, 20, 40, 60, 80, 100),
         las=1, lwd=0.8)
    title(main="Thresholding of mapping quality", cex.main=1.5)
    title(xlab="Quality threshold", line=2.5, cex.lab=1.2)
    title(ylab="Reads exceeding quality level, %", line=2.8, cex.lab=1.2)
    rect(xleft=0, ybottom=0, xright=10, ytop=100, border=NA,
         col=rgb(red=240, green=228, blue=66, maxColorValue=255, alpha=30))
    rect(xleft=10, ybottom=0, xright=20, ytop=100, border=NA,
         col=rgb(red=230, green=159, blue=0, maxColorValue=255, alpha=40))
    rect(xleft=20, ybottom=0, xright=30, ytop=100, border=NA,
         col=rgb(red=213, green=94, blue=0, maxColorValue=255, alpha=50))
    rect(xleft=30, ybottom=0, xright=40, ytop=100, border=NA,
         col=rgb(red=213, green=94, blue=0, maxColorValue=255, alpha=80))
    text(x=10, y=23, srt=90, cex=0.9,
         labels=paste(paste("MAPQ10",
                            round(x=length(x=x$MAPQs[x$MAPQs >= 10])/L * 100,
                                  digits=1),
                            sep=", "),
                      "%",
                      sep=""))
    text(x=20, y=23, srt=90, cex=0.9,
         labels=paste(paste("MAPQ20",
                            round(x=length(x=x$MAPQs[x$MAPQs >= 20])/L * 100,
                                  digits=1),
                            sep=", "),
                      "%",
                      sep=""))
    text(x=30, y=23, srt=90, cex=0.9,
         labels=paste(paste("MAPQ30",
                            round(x=length(x=x$MAPQs[x$MAPQs >= 30])/L * 100,
                                  digits=1),
                            sep=", "),
                      "%",
                      sep=""))
    ##  Plots of template lengths.
    if (x$TNRs > 1e6){
        p=sample(x=log10(x=x$TLENs + 1), size=1e6)
    }else{
        p=log10(x=x$TLENs + 1)
    }
    plot(0, type="n", axes=FALSE, ann=FALSE, bty="n",
         xlim=c(0, 2), ylim=c(0, ceiling(x=max(x=p))))
    vioplot(x=p,
            add=TRUE,
            col=rgb(red=86, green=180, blue=233, maxColorValue=255),
            rectCol=rgb(red=0, green=158, blue=115, maxColorValue=255),
            pchMed=0, colMed=NA,
            lwd=0.9,
            wex=1,
            frame.plot=FALSE)
    segments(x0=0.9, y0=median(x=p),
             x1=1.1, y1=median(x=p),
             col=rgb(red=240, green=228, blue=66, maxColorValue=255),
             lwd=3, lend=2)
    axis(side=1, at=c(0, 1, 2), labels=FALSE, lwd=0.8)
    axis(side=2,
         at=c(0.0, ceiling(x=max(x=p))/4 * 1:4),
         labels=round(x=c(0.0, ceiling(x=max(x=p))/4 * 1:4), digits=2),
         las=1, lwd=0.8)
    title(main="Distribution of template lengths", cex.main=1.5)
    title(xlab="Standard vioplot", line=1, cex.lab=1.2)
    title(ylab=expression("Length, log"[10] * " (nucleotides)"),
          line=2.8, cex.lab=1.2)
    p=log10(x=x$TLENs + 1)
    tr <- c(0.0, ceiling(x=max(x=p))/100 * 1:100)
    per <- data.frame(threshold=tr, percent=0)
    per[1, 2] <- 100
    L <- length(x=p)
    for (i in 1:(length(x=tr) - 1)){
        per[i + 1, 2] <- length(x=p[p >= tr[i]])/L * 100
    }
    plot(per, type="l", bty="n",
         xlim=c(0, ceiling(x=max(x=p))), ylim=c(0, 100),
         lwd=1,
         col=rgb(red=0, green=114, blue=178, maxColorValue=255),
         ann=FALSE, axes=FALSE)
    axis(side=1,
         at=c(0.0, ceiling(x=max(x=p))/4 * 1:4),
         labels=c(0.0, ceiling(x=max(x=p))/4 * 1:4),
         las=1, lwd=0.8)
    axis(side=2,
         at=c(0, 20, 40, 60, 80, 100),
         labels=c(0, 20, 40, 60, 80, 100),
         las=1, lwd=0.8)
    title(main="Thresholding of template lengths", cex.main=1.5)
    title(xlab=expression("Length, log"[10] * " (nucleotides)"),
          line=2.5, cex.lab=1.2)
    title(ylab="Templates exceeding length, %",
          line=2.8, cex.lab=1.2)
    rect(xleft=0, ybottom=0, xright=ceiling(x=max(x=p)), ytop=20, border=NA,
         col=rgb(red=240, green=228, blue=66, maxColorValue=255, alpha=30))
    rect(xleft=0, ybottom=20, xright=ceiling(x=max(x=p)), ytop=40, border=NA,
         col=rgb(red=230, green=159, blue=0, maxColorValue=255, alpha=40))
    rect(xleft=0, ybottom=40, xright=ceiling(x=max(x=p)), ytop=60, border=NA,
         col=rgb(red=213, green=94, blue=0, maxColorValue=255, alpha=50))
    rect(xleft=0, ybottom=60, xright=ceiling(x=max(x=p)), ytop=80, border=NA,
         col=rgb(red=213, green=94, blue=0, maxColorValue=255, alpha=70))
    rect(xleft=0, ybottom=80, xright=ceiling(x=max(x=p)), ytop=100, border=NA,
         col=rgb(red=213, green=94, blue=0, maxColorValue=255, alpha=90))
    ##  Title, subtitle, sub-subtitle & footnote.
    mtext(text="Mapping quality. Summary plots",
          side=3,
          line=3.5,
          cex=1.1,
          col="black",
          font=2,
          outer=TRUE)
    mtext(text=sprintf("Sample: %s", x$SN),
          side=3,
          line=2,
          cex=0.95,
          col="black",
          outer=TRUE)
    mtext(text=sprintf(paste("Date & Time:",
                             format(x=Sys.time(), "%Y-%m-%d %X"),
                             sep=" ")),
          side=3,
          line=0.5,
          cex=0.8,
          col="black",
          outer=TRUE)
    note <- "Generated by assessQAlignedReads() & reportQAlResults()"
    mtext(text=sprintf(note),
          side=1,
          line=2,
          cex=0.7,
          col="black",
          outer=TRUE)
    suppressMessages(expr=dev.off())
}
