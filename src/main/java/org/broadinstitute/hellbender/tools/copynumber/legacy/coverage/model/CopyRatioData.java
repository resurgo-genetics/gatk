package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.model;

import com.google.common.primitives.Doubles;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.Variance;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.TSVLocatableCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.DataCollection;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * {@link DataCollection} for the copy-ratio model containing the copy-ratio data grouped by segment.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class CopyRatioData implements DataCollection {
    private final CopyRatioCollection copyRatios;
    private final List<SimpleInterval> segments;
    private final double minLog2CopyRatioValue;
    private final double maxLog2CopyRatioValue;

    private final List<List<IndexedCopyRatio>> indexedCopyRatiosPerSegment = new ArrayList<>();

    CopyRatioData(final CopyRatioCollection copyRatios,
                  final List<SimpleInterval> segments) {
        this.copyRatios = Utils.nonNull(copyRatios);
        this.segments = Utils.nonEmpty(segments).stream().sorted(TSVLocatableCollection.LEXICOGRAPHICAL_ORDER_COMPARATOR).collect(Collectors.toList());

        final List<Double> log2CopyRatioValues = copyRatios.getLog2CopyRatioValues();
        minLog2CopyRatioValue = log2CopyRatioValues.stream().min(Double::compareTo).orElse(Double.NaN);
        maxLog2CopyRatioValue = log2CopyRatioValues.stream().max(Double::compareTo).orElse(Double.NaN);

        //construct list of lists of copy ratios with an index in order corresponding to that of segments;
        //segment assignment is based on midpoint of copy-ratio interval
        final OverlapDetector<CopyRatio> copyRatioMidpointOverlapDetector = copyRatios.getMidpointOverlapDetector();
        int index = 0;
        for (final SimpleInterval segment : segments) {
            final List<CopyRatio> copyRatiosInSegment = copyRatioMidpointOverlapDetector.getOverlaps(segment).stream()
                    .sorted(TSVLocatableCollection.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                    .collect(Collectors.toList());
            final int segmentStartIndex = index;
            final List<IndexedCopyRatio> indexedCopyRatiosInSegment = IntStream.range(0, copyRatiosInSegment.size()).boxed()
                    .map(i -> new IndexedCopyRatio(copyRatiosInSegment.get(i), segmentStartIndex + i))
                    .collect(Collectors.toList());
            indexedCopyRatiosPerSegment.add(indexedCopyRatiosInSegment);
            index += copyRatiosInSegment.size();
        }
    }

    int getNumSegments() {
        return segments.size();
    }

    int getNumPoints() {
        return copyRatios.size();
    }

    double getMinLog2CopyRatioValue() {
        return minLog2CopyRatioValue;
    }

    double getMaxLog2CopyRatioValue() {
        return maxLog2CopyRatioValue;
    }

    List<IndexedCopyRatio> getIndexedCopyRatiosInSegment(final int segment) {
        return Collections.unmodifiableList(indexedCopyRatiosPerSegment.get(segment));
    }

    //estimate global variance empirically by taking average of all per-segment variances
    double estimateVariance() {
        return IntStream.range(0, segments.size())
                .mapToDouble(s -> new Variance().evaluate(Doubles.toArray(
                        indexedCopyRatiosPerSegment.get(s).stream()
                                .map(IndexedCopyRatio::getLog2CopyRatioValue)
                                .collect(Collectors.toList()))))
                .filter(v -> !Double.isNaN(v))
                .average().orElse(Double.NaN);
    }

    //estimate segment means empirically by taking averages of log2 copy ratios in each segment
    CopyRatioState.SegmentMeans estimateSegmentMeans() {
        final List<Double> means = IntStream.range(0, segments.size()).boxed()
                .map(s -> new Mean().evaluate(Doubles.toArray(
                        indexedCopyRatiosPerSegment.get(s).stream()
                                .map(IndexedCopyRatio::getLog2CopyRatioValue)
                                .collect(Collectors.toList()))))
                .collect(Collectors.toList());
        return new CopyRatioState.SegmentMeans(means);
    }

    static final class IndexedCopyRatio extends CopyRatio {
        private final int index;

        private IndexedCopyRatio(final CopyRatio copyRatio,
                                 final int index) {
            super(copyRatio.getInterval(), copyRatio.getLog2CopyRatioValue());
            this.index = index;
        }

        int getIndex() {
            return index;
        }
    }
}
