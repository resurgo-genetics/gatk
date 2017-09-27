#NOTE: the Java wrapper for this script first sources CNVPlottingLibrary.R
options(error = quote({dump.frames(dumpto = "plotting_dump", to.file = TRUE); q(status = 1)}))    # Useful for debugging

library(optparse)
option_list = list(
    make_option(c("--sample_name", "-sample_name"), dest="sample_name", action="store"),
    make_option(c("--standardized_file", "-standardized_file"), dest="standardized_file", action="store"),
    make_option(c("--denoised_file", "-denoised_file"), dest="denoised_file", action="store"),
    make_option(c("--contig_names", "-contig_names"), dest="contig_names", action="store"),         #string with elements separated by "CONTIG_DELIMITER"
    make_option(c("--contig_lengths", "-contig_lengths"), dest="contig_lengths", action="store"),   #string with elements separated by "CONTIG_DELIMITER"
    make_option(c("--output_dir", "-output_dir"), dest="output_dir", action="store"),
    make_option(c("--output_prefix", "-output_prefix"), dest="output_prefix", action="store"))

opt = parse_args(OptionParser(option_list=option_list))

sample_name = opt[["sample_name"]]
standardized_file = opt[["standardized_file"]]
denoised_file = opt[["denoised_file"]]
contig_names_string = opt[["contig_names"]]
contig_lengths_string = opt[["contig_lengths"]]
output_dir = opt[["output_dir"]]
output_prefix = opt[["output_prefix"]]

#check that input files exist; if not, quit with error code that GATK will pick up
if (!all(file.exists(c(standardized_file, denoised_file)))) {
    quit(save="no", status=1, runLast=FALSE)
}

contig_names = as.list(strsplit(contig_names_string, "CONTIG_DELIMITER")[[1]])
contig_lengths = as.list(strsplit(contig_lengths_string, "CONTIG_DELIMITER")[[1]])
contig_ends = cumsum(contig_lengths)
contig_starts = c(0, head(contig_ends, -1))

CalculateQc = function(dat) {
    return(median(abs(diff(dat))))
}

#plotting is extracted to a function for debugging purposes
plot_denoised_copy_ratios = function(standardized_file, denoised_file, contig_names, output_dir, output_prefix) {
    #read in files and extract needed data
    standardized = read.table(standardized_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
    denoised = read.table(denoised_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)

    #transform to linear copy ratio
    standardized$COPY_RATIO = 2^standardized$LOG2_COPY_RATIO
    denoised$COPY_RATIO = 2^denoised$LOG2_COPY_RATIO

    #write the QC files
    preQc = CalculateQc(standardized$COPY_RATIO)
    postQc = CalculateQc(denoised$COPY_RATIO)
    write.table(round(preQc, 3), file.path(output_dir, paste(output_prefix, "_preQc.txt", sep="")), col.names=FALSE, row.names=FALSE)
    write.table(round(postQc, 3), file.path(output_dir, paste(output_prefix, "_postQc.txt", sep="")), col.names=FALSE, row.names=FALSE)
    write.table(round(preQc - postQc, 3), file.path(output_dir, paste(output_prefix, "_dQc.txt", sep="")), col.names=FALSE, row.names=FALSE)
    write.table(round((preQc - postQc) / preQc, 3), file.path(output_dir, paste(output_prefix, "_scaled_dQc.txt", sep="")), col.names=FALSE, row.names=FALSE)

    #plot standardized and denoised copy ratio on top of each other
    pre_color_blue="#3B5DFF"
    post_color_green="#4FC601"
    #plot over full range
    plot_before_after_full_file = file.path(output_dir, paste(output_prefix, "_Before_After.png", sep=""))
    png(plot_before_after_full_file, 12, 7, units="in", type="cairo", res=300, bg="white")
    par(mfrow=c(2,1), cex=0.75, las=1)
    SetUpPlot("Standardized copy ratio", 0, max(standardized$COPY_RATIO), paste("After standardization, QC = ", round(preQc, 3), sep=""), contig_names, contig_starts, contig_ends, FALSE)
    PlotCopyRatio(standardized, pre_color_blue, contig_names, contig_starts)
    SetUpPlot("Denoised copy ratio", 0, max(denoised$COPY_RATIO), paste("After denoising, QC = ", round(postQc, 3), sep=""), contig_names, contig_starts, contig_ends, TRUE)
    PlotCopyRatio(denoised, post_color_green, contig_names, contig_starts)
    dev.off()
    #plot up to CR = 4
    plot_before_after_CR_lim_file = file.path(output_dir, paste(output_prefix, "_Before_After_CR_Lim_4.png", sep=""))
    png(plot_before_after_CR_lim_file, 12, 7, units="in", type="cairo", res=300, bg="white")
    par(mfrow=c(2,1), cex=0.75, las=1)
    SetUpPlot("Standardized copy ratio", 0, 4, paste("After standardization, QC = ", round(preQc, 3), sep=""), contig_names, contig_starts, contig_ends, FALSE)
    PlotCopyRatio(standardized, pre_color_blue, contig_names, contig_starts)
    SetUpPlot("Denoised copy ratio", 0, 4, paste("After denoising, QC = ", round(postQc, 3), sep=""), contig_names, contig_starts, contig_ends, TRUE)
    PlotCopyRatio(denoised, post_color_green, contig_names, contig_starts)
    dev.off()

    #check for created files and quit with error code if not found
    if (!all(file.exists(c(plot_before_after_full_file, plot_before_after_CR_lim_file)))) {
        quit(save="no", status=1, runLast=FALSE)
    }
}

plot_denoised_copy_ratios(standardized_file, denoised_file, contig_names, output_dir, output_prefix)