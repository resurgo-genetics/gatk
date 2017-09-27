## Contig Background for data
SetUpPlot = function(y.lab, y.min, y.max, x.lab, contig_names, contig_starts, contig_ends, do_label_contigs) {
    num_contigs = length(contig_names)
    contig_centers = (contig_starts + contig_ends) / 2
    genome_length = contig_ends[num_contigs]
    suppressWarnings(par(mar=c(3.1, 3.6, 0.1, 0), mgp=c(2, -0.2, -1.1)))
    plot(0, type="n", bty="n", xlim=c(0, genome_length), ylim=c(0, y.max), xlab="", ylab="", main="", xaxt="n")
    mtext(side=1, line=1.5, x.lab, cex=0.75, outer=FALSE)
    mtext(side=2.2, line=1.5, y.lab, cex=0.75, las=FALSE, outer=FALSE)

    if (do_label_contigs) {
        mtext(text=contig_names[1:num_contigs], side=1, line=ifelse(c(1:num_contigs) %% 2 == 1, -0.45, 0.0),
        at = contig_centers[1:num_contigs], las=1, cex=par("cex.axis") * par("cex") * 0.7)
    }

    for (i in 1:num_contigs) {
        use.col = ifelse(i %% 2 == 1, "grey90", "white")
        rect(xleft=contig_starts[i], ybottom=y.min, xright=contig_ends[i], ytop=y.max, col=use.col, border=NA)
    }
}

PlotCopyRatio = function(copy_ratio_df, color, contig_names, contig_starts) {
    genomic_coordinates = contig_starts[match(copy_ratio_df$CONTIG, contig_names)] + copy_ratio_df$END
    points(x=genomic_coordinates, y=copy_ratio_df$COPY_RATIO, col=color, pch=16, cex=0.2)
}