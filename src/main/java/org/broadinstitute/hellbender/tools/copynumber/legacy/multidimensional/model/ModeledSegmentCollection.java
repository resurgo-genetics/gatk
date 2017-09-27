package org.broadinstitute.hellbender.tools.copynumber.legacy.multidimensional.model;

import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.TSVLocatableCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class ModeledSegmentCollection extends TSVLocatableCollection<ModeledSegment> {
    public static final String DOUBLE_FORMAT = "%6.6f";

    enum ModeledSegmentTableColumn {
        CONTIG, 
        START, 
        END,
        NUM_POINTS_COPY_RATIO,
        NUM_POINTS_ALLELE_FRACTION,
        LOG2_COPY_RATIO_POSTERIOR_MODE,
        LOG2_COPY_RATIO_POSTERIOR_10,
        LOG2_COPY_RATIO_POSTERIOR_50,
        LOG2_COPY_RATIO_POSTERIOR_90,
        MINOR_ALLELE_FRACTION_POSTERIOR_MODE,
        MINOR_ALLELE_FRACTION_POSTERIOR_10,
        MINOR_ALLELE_FRACTION_POSTERIOR_50,
        MINOR_ALLELE_FRACTION_POSTERIOR_90;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static final Function<DataLine, ModeledSegment> MODELED_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION = dataLine -> {
        final String contig = dataLine.get(ModeledSegmentTableColumn.CONTIG);
        final int start = dataLine.getInt(ModeledSegmentTableColumn.START);
        final int end = dataLine.getInt(ModeledSegmentTableColumn.END);
        final int numPointsCopyRatio = dataLine.getInt(ModeledSegmentTableColumn.NUM_POINTS_COPY_RATIO);
        final int numPointsAlleleFraction = dataLine.getInt(ModeledSegmentTableColumn.NUM_POINTS_ALLELE_FRACTION);
        final double log2CopyRatioPosteriorMode = dataLine.getDouble(ModeledSegmentTableColumn.LOG2_COPY_RATIO_POSTERIOR_MODE);
        final double log2CopyRatioPosterior10 = dataLine.getDouble(ModeledSegmentTableColumn.LOG2_COPY_RATIO_POSTERIOR_10);
        final double log2CopyRatioPosterior50 = dataLine.getDouble(ModeledSegmentTableColumn.LOG2_COPY_RATIO_POSTERIOR_50);
        final double log2CopyRatioPosterior90 = dataLine.getDouble(ModeledSegmentTableColumn.LOG2_COPY_RATIO_POSTERIOR_90);
        final double minorAlleleFractionPosteriorMode = dataLine.getDouble(ModeledSegmentTableColumn.MINOR_ALLELE_FRACTION_POSTERIOR_MODE);
        final double minorAlleleFractionPosterior10 = dataLine.getDouble(ModeledSegmentTableColumn.MINOR_ALLELE_FRACTION_POSTERIOR_10);
        final double minorAlleleFractionPosterior50 = dataLine.getDouble(ModeledSegmentTableColumn.MINOR_ALLELE_FRACTION_POSTERIOR_50);
        final double minorAlleleFractionPosterior90 = dataLine.getDouble(ModeledSegmentTableColumn.MINOR_ALLELE_FRACTION_POSTERIOR_90);
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        return new ModeledSegment(interval, numPointsCopyRatio, numPointsAlleleFraction, 
                log2CopyRatioPosteriorMode, log2CopyRatioPosterior10, log2CopyRatioPosterior50, log2CopyRatioPosterior90,
                minorAlleleFractionPosteriorMode, minorAlleleFractionPosterior10, minorAlleleFractionPosterior50, minorAlleleFractionPosterior90);
    };

    private static final BiConsumer<ModeledSegment, DataLine> MODELED_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER = (alleleFractionSegment, dataLine) ->
            dataLine.append(alleleFractionSegment.getContig())
                    .append(alleleFractionSegment.getStart())
                    .append(alleleFractionSegment.getEnd())
                    .append(alleleFractionSegment.getNumPointsCopyRatio())
                    .append(alleleFractionSegment.getNumPointsAlleleFraction())
                    .append(String.format(DOUBLE_FORMAT, alleleFractionSegment.getLog2CopyRatioPosteriorMode()))
                    .append(String.format(DOUBLE_FORMAT, alleleFractionSegment.getLog2CopyRatioPosterior10()))
                    .append(String.format(DOUBLE_FORMAT, alleleFractionSegment.getLog2CopyRatioPosterior50()))
                    .append(String.format(DOUBLE_FORMAT, alleleFractionSegment.getLog2CopyRatioPosterior90()))
                    .append(String.format(DOUBLE_FORMAT, alleleFractionSegment.getMinorAlleleFractionPosteriorMode()))
                    .append(String.format(DOUBLE_FORMAT, alleleFractionSegment.getMinorAlleleFractionPosterior10()))
                    .append(String.format(DOUBLE_FORMAT, alleleFractionSegment.getMinorAlleleFractionPosterior50()))
                    .append(String.format(DOUBLE_FORMAT, alleleFractionSegment.getMinorAlleleFractionPosterior90()));

    public ModeledSegmentCollection(final File inputFile) {
        super(inputFile, ModeledSegmentTableColumn.COLUMNS, MODELED_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION, MODELED_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }

    public ModeledSegmentCollection(final String sampleName,
                                    final List<ModeledSegment> modeledSegments) {
        super(sampleName, modeledSegments, ModeledSegmentTableColumn.COLUMNS, MODELED_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION, MODELED_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }
}