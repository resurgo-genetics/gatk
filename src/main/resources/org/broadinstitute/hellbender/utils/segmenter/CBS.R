# Useful for debugging
#options(error = quote({dump.frames(dumpto = "CBS_dump", to.file = TRUE); q()}))

# Library used for segmentation
library(DNAcopy)
library(naturalsort)

library(optparse)
option_list <- list(
    make_option(c("--sample_name", "-sample_name"), dest="sample_name", action="store"),
    make_option(c("--targets_file", "-targets_file"), dest="targets_file", action="store"),
    make_option(c("--output_file", "-output_file"), dest="output_file", action="store"),
    make_option(c("--log2_input", "-log"), dest="log2_input", action="store"),
    make_option(c("--min_width", "-mw"), dest="min_width", action="store", default=2),
    make_option(c("--weights_file", "-w"), dest="weights_file", action="store", default=NULL),
    make_option(c("--alpha" , "-alpha" ), dest="alpha" , action="store", type="double"),
    make_option(c("--nperm" , "-nperm" ), dest="nperm" , action="store", type="integer"),
    make_option(c("--pmethod" , "-pmethod" ), dest="pmethod" , action="store"),
    make_option(c("--kmax" , "-kmax" ), dest="kmax" , action="store", type="integer"),
    make_option(c("--nmin" , "-nmin" ), dest="nmin" , action="store", type="integer"),
    make_option(c("--eta" , "-eta" ), dest="eta" , action="store", type="double"),
    make_option(c("--trim" , "-trim" ), dest="trim" , action="store", type="double"),
    make_option(c("--undosplits" , "-undosplits" ), dest="undosplits" , action="store"),
    make_option(c("--undoprune" , "-undoprune" ), dest="undoprune" , action="store"),
    make_option(c("--undoSD" , "-undoSD" ), dest="undoSD" , action="store", type="integer")
    )

opt <- parse_args(OptionParser(option_list=option_list))
print(opt)

sample_name=opt[["sample_name"]]
tn_file=opt[["targets_file"]]
output_file=opt[["output_file"]]
log_input=as.logical(opt[["log2_input"]])
min_width=opt[["min_width"]]
weights_file = opt[["weights_file"]]
alpha=opt[["alpha" ]]
nperm=opt[["nperm" ]]
pmethod=opt[["pmethod" ]]
kmax=opt[["kmax" ]]
nmin=opt[["nmin" ]]
eta=opt[["eta" ]]
trim=opt[["trim" ]]
undosplits=opt[["undosplits" ]]
undoprune=opt[["undoprune" ]]
undoSD=opt[["undoSD" ]]

# segment(x, weights = NULL, alpha = 0.01, nperm = 10000, p.method =
#  c("hybrid", "perm"), min.width=2, kmax=25, nmin=200,
#  eta=0.05, sbdry=NULL, trim = 0.025, undo.splits =
#  c("none", "prune", "sdundo"), undo.prune=0.05,
#  undo.SD=3, verbose=1)

# Use a function for debugging purposes
segment_data = function(sample_name, tn_file, output_file, log_input, weights_file, min_width, alpha, nperm, pmethod, kmax, nmin, eta, trim, undosplits, undoprune, undoSD) {
	# Read in file and extract needed data
	tn = read.table(tn_file, sep="\t", stringsAsFactors=FALSE, header=TRUE, check.names=FALSE)
	contig = tn[,"contig"]
	pos = tn[,"stop"]
	# Ability to switch between copy-ratio and log2 copy-ratio
	if (log_input) {
	    dat = tn[,sample_name]
	} else {
	    dat = log2(tn[,sample_name])
	}

	if (is.null(weights_file)) {
		weights = NULL
	} else {
		# weights needs to be in a row vector, so transpose (t(...)) is needed
		weights = t(read.table(weights_file, sep="\t", header=FALSE))
	}

	# Create CNA object
	cna_dat = CNA(dat, contig, pos, data.type="logratio", sampleid=sample_name)

	# Perform segmentation
	set.seed(25)

	# segment has an issue with passing in weights of NULL.  Better to not specify.
	if (!is.null(weights_file)) {
		segmented = segment(smooth.CNA(cna_dat), min.width=min_width, weights=weights,alpha=alpha, nperm=nperm, p.method=pmethod, kmax=kmax, nmin=nmin, eta=eta, trim=trim, undo.splits=undosplits, undo.prune=undoprune, undo.SD=undoSD)$output
	} else {
		segmented = segment(smooth.CNA(cna_dat), min.width=min_width, alpha=alpha, nperm=nperm, p.method=pmethod, kmax=kmax, nmin=nmin, eta=eta, trim=trim, undo.splits=undosplits, undo.prune=undoprune, undo.SD=undoSD)$output
	}

	# Ensure that there are no too-small values which will be problematic for downstream tools.
	segmented[,"seg.mean"] = 2^segmented[,"seg.mean"]
	segmented[segmented[,"seg.mean"]<.Machine$double.eps,"seg.mean"] = .Machine$double.eps

	# Convention for column names
	colnames(segmented) = c("Sample", "Chromosome", "Start", "End", "Num_Probes", "Segment_Mean")

	# Undo conversion of sample_name to R format (retain dashes, spaces, etc.)
	segmented[,"Sample"] = sample_name

	# Change segment start to start of first target in segment (rather than end of first target in segment)
	segmented$rowID = 1:nrow(segmented)
	first_targets = merge(segmented[,c("rowID","Chromosome","Start")], tn[,c("contig","start","stop")], by.x=c("Chromosome","Start"), by.y=c("contig","stop"))
	# merge does not maintain order; restore order of original segmented dataframe
	first_targets = first_targets[order(first_targets$rowID), ]
	segmented[,"Start"] = first_targets[,"start"]
	segmented$rowID = NULL

    # Order based on contig (already ordered based on start)
    sorting = unique(naturalsort(segmented[,"Chromosome"]))
    segmented$Chromosome=factor(segmented$Chromosome, levels=sorting)
    segmented = segmented[order(segmented[,"Chromosome"]),]

	# Output seg file
	print(paste("Writing segment file: ", output_file, sep=""))
	write.table(segmented, file=output_file, sep="\t", quote=FALSE, row.names=FALSE)
}

segment_data(sample_name, tn_file, output_file, log_input, weights_file, min_width, alpha, nperm, pmethod, kmax, nmin, eta, trim, undosplits, undoprune, undoSD)
