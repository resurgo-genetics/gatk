package org.broadinstitute.hellbender.tools.copynumber.legacy.multidimensional.segmentation;

import com.google.common.primitives.Doubles;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.segmentation.AlleleFractionSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.TSVLocatableCollection;
import org.broadinstitute.hellbender.tools.exome.*;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Represents a legacy allele-fraction segmentation.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CRAFSegmentCollection extends TSVLocatableCollection<CRAFSegment> {
    enum CRAFSegmentTableColumn {
        CONTIG,
        START,
        END,
        NUM_POINTS_ALLELE_FRACTION,
        NUM_POINTS_COPY_RATIO,
        MEAN_MINOR_ALLELE_FRACTION,
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
        final double meanMinorAlleleFraction = dataLine.getDouble(CRAFSegmentTableColumn.MEAN_MINOR_ALLELE_FRACTION);
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        return new CRAFSegment(interval, numPointsCopyRatio, numPointsAlleleFraction, meanLog2CopyRatio, meanMinorAlleleFraction);
    };

    private static final BiConsumer<CRAFSegment, DataLine> CRAF_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER = (alleleFractionSegment, dataLine) ->
            dataLine.append(alleleFractionSegment.getContig())
                    .append(alleleFractionSegment.getStart())
                    .append(alleleFractionSegment.getEnd())
                    .append(alleleFractionSegment.getNumPointsCopyRatio())
                    .append(alleleFractionSegment.getNumPointsAlleleFraction())
                    .append(alleleFractionSegment.getMeanLog2CopyRatio())
                    .append(alleleFractionSegment.getMeanMinorAlleleFraction());

    public CRAFSegmentCollection(final File inputFile) {
        super(inputFile, CRAFSegmentTableColumn.COLUMNS, CRAF_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION, CRAF_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }

    public CRAFSegmentCollection(final String sampleName,
                                 final List<CRAFSegment> crafSegments) {
        super(sampleName, crafSegments, CRAFSegmentTableColumn.COLUMNS, CRAF_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION, CRAF_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }

    public static CRAFSegmentCollection unionAndMergeSmallSegments(final CopyRatioSegmentCollection copyRatioSegments,
                                                                   final CopyRatioCollection denoisedCopyRatios,
                                                                   final AlleleFractionSegmentCollection alleleFractionSegments,
                                                                   final AllelicCountCollection allelicCounts,
                                                                   final int numCopyRatioIntervalsSmallSegmentThreshold) {
        Utils.validateArg(!(copyRatioSegments == null && alleleFractionSegments == null),
                "Must provide at least a copy-ratio segmentation or an allele-fraction segmentation.");
        ParamUtils.isPositiveOrZero(numCopyRatioIntervalsSmallSegmentThreshold,
                "Threshold number of copy-ratio intervals for small-segment merging must be non-negative.");

        final String sampleName;
        final List<CRAFSegment> crafSegments;

        if (alleleFractionSegments == null) {
            Utils.nonNull(denoisedCopyRatios);
            sampleName = copyRatioSegments.getSampleName();
            crafSegments = copyRatioSegments.getRecords().stream()
                    .map(s -> new CRAFSegment(s.getInterval(), s.getNumPoints(), 0, s.getMeanLog2CopyRatio(), Double.NaN))
                    .collect(Collectors.toList());
        } else if (copyRatioSegments == null) {
            Utils.nonNull(allelicCounts);
            sampleName = alleleFractionSegments.getSampleName();
            crafSegments = alleleFractionSegments.getRecords().stream()
                    .map(s -> new CRAFSegment(s.getInterval(), 0, s.getNumPoints(), Double.NaN, s.getMeanMinorAlleleFraction()))
                    .collect(Collectors.toList());
        } else {
            Utils.validateArg(copyRatioSegments.getSampleName().equals(alleleFractionSegments.getSampleName()),
                    "Sample names from copy-ratio segmentation and allele-fraction segmentation must match.");
            sampleName = copyRatioSegments.getSampleName();
            //use old code for segment union
            final Genome genome = new Genome(
                    new ReadCountCollection(
                            denoisedCopyRatios.getIntervals().stream().map(Target::new).collect(Collectors.toList()),
                            Collections.singletonList(sampleName),
                            new Array2DRowRealMatrix(Doubles.toArray(denoisedCopyRatios.getCopyRatioValues()))),
                    allelicCounts.getRecords().stream()
                            .map(ac -> new org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount(ac.getInterval(), ac.getRefReadCount(), ac.getAltReadCount()))
                            .collect(Collectors.toList()));
            final List<SimpleInterval> unionedSegments = SegmentUtils.unionSegments(copyRatioSegments.getIntervals(), alleleFractionSegments.getIntervals(), genome);
            final SegmentedGenome segmentedGenomeWithSmallSegments = new SegmentedGenome(unionedSegments, genome);
            final SegmentedGenome segmentedGenome = segmentedGenomeWithSmallSegments.mergeSmallSegments(numCopyRatioIntervalsSmallSegmentThreshold);
            final OverlapDetector<CopyRatio> copyRatioOverlapDetector = OverlapDetector.create(denoisedCopyRatios.getRecords());
            final OverlapDetector<AllelicCount> allelicCountOverlapDetector = OverlapDetector.create(allelicCounts.getRecords());
            crafSegments = segmentedGenome.getSegments().stream()
                    .map(s -> new CRAFSegment(
                            s,
                            copyRatioOverlapDetector.getOverlaps(s).stream()
                                    .map(CopyRatio::getLog2CopyRatioValue)
                                    .collect(Collectors.toList()),
                            new ArrayList<>(allelicCountOverlapDetector.getOverlaps(s))))
                    .collect(Collectors.toList());
        }

        return new CRAFSegmentCollection(sampleName, crafSegments);
    }
}