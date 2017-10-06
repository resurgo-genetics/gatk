package org.broadinstitute.hellbender.tools.copynumber.multidimensional.segmentation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.allelic.segmentation.AlleleFractionSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.coverage.segmentation.CopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.LocatableCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Stream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CRAFSegmentCollection extends LocatableCollection<CRAFSegment> {
    private static final Logger logger = LogManager.getLogger(CRAFSegmentCollection.class);

    enum CRAFSegmentTableColumn {
        CONTIG,
        START,
        END,
        NUM_POINTS_COPY_RATIO,
        NUM_POINTS_ALLELE_FRACTION,
        MEAN_LOG2_COPY_RATIO;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static final Function<DataLine, CRAFSegment> CRAF_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION = dataLine -> {
        final String contig = dataLine.get(CRAFSegmentTableColumn.CONTIG);
        final int start = dataLine.getInt(CRAFSegmentTableColumn.START);
        final int end = dataLine.getInt(CRAFSegmentTableColumn.END);
        final int numPointsCopyRatio = dataLine.getInt(CRAFSegmentTableColumn.NUM_POINTS_COPY_RATIO);
        final int numPointsAlleleFraction = dataLine.getInt(CRAFSegmentTableColumn.NUM_POINTS_ALLELE_FRACTION);
        final double meanLog2CopyRatio = dataLine.getDouble(CRAFSegmentTableColumn.MEAN_LOG2_COPY_RATIO);
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        return new CRAFSegment(interval, numPointsCopyRatio, numPointsAlleleFraction, meanLog2CopyRatio);
    };

    private static final BiConsumer<CRAFSegment, DataLine> CRAF_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER = (alleleFractionSegment, dataLine) ->
            dataLine.append(alleleFractionSegment.getContig())
                    .append(alleleFractionSegment.getStart())
                    .append(alleleFractionSegment.getEnd())
                    .append(alleleFractionSegment.getNumPointsCopyRatio())
                    .append(alleleFractionSegment.getNumPointsAlleleFraction())
                    .append(alleleFractionSegment.getMeanLog2CopyRatio());

    public CRAFSegmentCollection(final File inputFile) {
        super(inputFile, CRAFSegmentTableColumn.COLUMNS, CRAF_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION, CRAF_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }

    public CRAFSegmentCollection(final String sampleName,
                                 final List<CRAFSegment> crafSegments) {
        super(sampleName, crafSegments, CRAFSegmentTableColumn.COLUMNS, CRAF_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION, CRAF_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }

    public static CRAFSegmentCollection unionSegments(final CopyRatioSegmentCollection copyRatioSegments,
                                                      final CopyRatioCollection denoisedCopyRatios,
                                                      final AlleleFractionSegmentCollection alleleFractionSegments,
                                                      final AllelicCountCollection allelicCounts) {
        final String sampleName;
        final List<CRAFSegment> crafSegments;

        Utils.validateArg(Stream.of(
                Utils.nonNull(copyRatioSegments).getSampleName(),
                Utils.nonNull(denoisedCopyRatios).getSampleName(),
                Utils.nonNull(alleleFractionSegments).getSampleName(),
                Utils.nonNull(allelicCounts).getSampleName())
                        .distinct().count() == 1,
                "Sample names from all inputs must match.");
        sampleName = copyRatioSegments.getSampleName();

        //union copy-ratio and allele-fraction segments
        logger.info(String.format("Combining %d copy-ratio segments and %d allele-fraction segments...",
                copyRatioSegments.getRecords().size(), alleleFractionSegments.getRecords().size()));
        crafSegments = new SegmentUnioner(copyRatioSegments, denoisedCopyRatios, alleleFractionSegments, allelicCounts)
                .constructUnionedCRAFSegments();
        logger.info(String.format("After combining segments, %d segments remain...", crafSegments.size()));

        return new CRAFSegmentCollection(sampleName, crafSegments);
    }
}