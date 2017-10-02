package org.broadinstitute.hellbender.tools.copynumber.legacy.multidimensional.segmentation;

import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.segmentation.AlleleFractionSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.TSVLocatableCollection;
import org.broadinstitute.hellbender.tools.exome.Genome;
import org.broadinstitute.hellbender.tools.exome.SegmentedGenome;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.Arrays;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CRAFSegmentCollection extends TSVLocatableCollection<CRAFSegment> {
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
        Utils.validateArg(!(copyRatioSegments == null && alleleFractionSegments == null),
                "Must provide at least a copy-ratio segmentation or an allele-fraction segmentation.");

        final String sampleName;
        final List<CRAFSegment> crafSegments;

        if (alleleFractionSegments == null) {
            Utils.nonNull(denoisedCopyRatios);
            Utils.validateArg(copyRatioSegments.getSampleName().equals(denoisedCopyRatios.getSampleName()),
                    "Sample names from copy-ratio segmentation and denoised copy ratios must match.");
            sampleName = copyRatioSegments.getSampleName();
            crafSegments = copyRatioSegments.getRecords().stream()
                    .map(s -> new CRAFSegment(s.getInterval(), s.getNumPoints(), 0, s.getMeanLog2CopyRatio()))
                    .collect(Collectors.toList());
        } else if (copyRatioSegments == null) {
            Utils.nonNull(allelicCounts);
            Utils.validateArg(alleleFractionSegments.getSampleName().equals(allelicCounts.getSampleName()),
                    "Sample names from allele-fraction segmentation and allelic counts must match.");
            sampleName = alleleFractionSegments.getSampleName();
            crafSegments = alleleFractionSegments.getRecords().stream()
                    .map(s -> new CRAFSegment(s.getInterval(), 0, s.getNumPoints(), Double.NaN))
                    .collect(Collectors.toList());
        } else {
            Utils.validateArg(Stream.of(
                    copyRatioSegments.getSampleName(),
                    denoisedCopyRatios.getSampleName(),
                    alleleFractionSegments.getSampleName(),
                    allelicCounts.getSampleName())
                    .distinct().count() == 1,
                    "Sample names from all inputs must match.");
            sampleName = copyRatioSegments.getSampleName();

            //union copy-ratio and allele-fraction segments
            logger.info(String.format("Combining %d copy-ratio segments and %d allele-fraction segments...",
                    copyRatioSegments.getRecords().size(), alleleFractionSegments.getRecords().size()));
            crafSegments = new SegmentUnioner(copyRatioSegments, denoisedCopyRatios, alleleFractionSegments, allelicCounts)
                    .constructUnionedCRAFSegments();
            logger.info(String.format("After combining segments, %d segments remain...", crafSegments.size()));
        }
        return new CRAFSegmentCollection(sampleName, crafSegments);
    }

    public SegmentedGenome convertToSegmentedGenome(final CopyRatioCollection denoisedCopyRatios,
                                                    final AllelicCountCollection allelicCounts) {
        final Genome genome = new Genome(denoisedCopyRatios, allelicCounts);
        return new SegmentedGenome(getIntervals(), genome);
    }
}