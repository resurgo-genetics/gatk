#NOTE: the Java wrapper for this script first sources CNVPlottingLibrary.R
options(error = quote({dump.frames(dumpto = "plotting_dump", to.file = TRUE); q(status = 1)}))    # Useful for debugging

library(optparse)
option_list = list(
make_option(c("--sample_name", "-sample_name"), dest="sample_name", action="store"),
make_option(c("--denoised_file", "-denoised_file"), dest="denoised_file", action="store"),
make_option(c("--allelic_counts_file", "-allelic_counts_file"), dest="allelic_counts_file", action="store"),
make_option(c("--modeled_segments_file", "-modeled_segments_file"), dest="modeled_segments_file", action="store"),
make_option(c("--contig_names", "-contig_names"), dest="contig_names", action="store"),         #string with elements separated by "CONTIG_DELIMITER"
make_option(c("--contig_lengths", "-contig_lengths"), dest="contig_lengths", action="store"),   #string with elements separated by "CONTIG_DELIMITER"
make_option(c("--output_dir", "-output_dir"), dest="output_dir", action="store"),
make_option(c("--output_prefix", "-output_prefix"), dest="output_prefix", action="store"))
opt = parse_args(OptionParser(option_list=option_list))

sample_name = opt[["sample_name"]]
denoised_file = opt[["denoised_file"]]
allelic_counts_file = opt[["allelic_counts_file"]]
modeled_segments_file = opt[["modeled_segments_file"]]
contig_names_string = opt[["contig_names"]]
contig_lengths_string = opt[["contig_lengths"]]
output_dir = opt[["output_dir"]]
output_prefix = opt[["output_prefix"]]

#check that input files exist; if not, quit with error code that GATK will pick up
if (!all(file.exists(c(allelic_counts_file, denoised_file, modeled_segments_file)))) {
    quit(save="no", status=1, runLast=FALSE)
}

contig_names = as.list(strsplit(contig_names_string, "CONTIG_DELIMITER")[[1]])
contig_lengths = as.list(strsplit(contig_lengths_string, "CONTIG_DELIMITER")[[1]])
contig_ends = cumsum(contig_lengths)
contig_starts = c(0, head(contig_ends, -1))

#plotting is extracted to a function for debugging purposes
write_modeled_segments_plot = function(sample_name, allelic_counts_file, denoised_file, modeled_segments_file, contig_names, contig_lengths, output_dir, output_prefix) {
    #set up copy ratio, allelics, and modeled_segments data frames
    copy_ratio = read.table(denoised_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
    allelic_counts = read.table(allelic_counts_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
    modeled_segments = read.table(modeled_segments_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)

    #transform to linear copy ratio
    copy_ratio$COPY_RATIO = 2^copy_ratio$LOG2_COPY_RATIO

    #plot CR and AAF data and segment posteriors
    plot_file = file.path(output_dir, paste(output_prefix, "_modeled_segments.png", sep=""))
    png(plot_file, 12, 7, units="in", type="cairo", res=300, bg="white")
    par(mfrow=c(2,1), cex=0.75, las=1)
    SetUpPlot("Denoised copy ratio", 0, 4, "Contig", contig_names, contig_starts, contig_ends, TRUE)
    PlotCopyRatioWithModeledSegments(copy_ratio, modeled_segments, contig_names, contig_starts, TRUE)
    SetUpPlot("Alternate-allele fraction", 0, 0.5, "Contig", contig_names, contig_starts, contig_ends, TRUE)
    PlotAlternateAlleleFractionWithModeledSegments(allelic_counts, modeled_segments, contig_names, contig_starts)
    dev.off()

    #check for created file and quit with error code if not found
    if (!file.exists(plot_file)) {
        quit(save="no", status=1, runLast=FALSE)
    }
}

write_modeled_segments_plot(sample_name, allelic_counts_file, denoised_file, modeled_segments_file, contig_names, contig_lengths, output_dir, output_prefix)