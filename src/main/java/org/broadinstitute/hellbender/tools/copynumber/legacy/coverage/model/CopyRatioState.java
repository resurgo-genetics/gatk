package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.model;

import org.broadinstitute.hellbender.utils.mcmc.Parameter;
import org.broadinstitute.hellbender.utils.mcmc.ParameterizedState;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;

/**
 * The state of the copy-ratio model, containing: <p>
 *      1.  the global variance <p>
 *      2.  the global outlier probability <p>
 *      3.  log2 mean copy ratios for each segment <p>
 *      4.  outlier indicators for each copy-ratio interval <p>
 * <p>
 * See docs/CNVs/CNV-methods.pdf for details.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class CopyRatioState extends ParameterizedState<CopyRatioParameter> {
    static final class SegmentMeans extends ArrayList<Double> {
        private static final long serialVersionUID = 951753L;

        SegmentMeans(final List<Double> segmentMeans) {
            super(new ArrayList<>(segmentMeans));
        }
    }

    static final class OutlierIndicators extends ArrayList<Boolean> {
        private static final long serialVersionUID = 357159L;

        OutlierIndicators(final List<Boolean> outlierIndicators) {
            super(new ArrayList<>(outlierIndicators));
        }
    }

    CopyRatioState(final double variance,
                   final double outlierProbability,
                   final SegmentMeans segmentMeans,
                   final OutlierIndicators outlierIndicators) {
        super(Arrays.asList(
                new Parameter<>(CopyRatioParameter.VARIANCE, variance),
                new Parameter<>(CopyRatioParameter.OUTLIER_PROBABILITY, outlierProbability),
                new Parameter<>(CopyRatioParameter.SEGMENT_MEANS, segmentMeans),
                new Parameter<>(CopyRatioParameter.OUTLIER_INDICATORS, outlierIndicators)));
    }

    double variance() {
        return get(CopyRatioParameter.VARIANCE, Double.class);
    }

    double outlierProbability() {
        return get(CopyRatioParameter.OUTLIER_PROBABILITY, Double.class);
    }

    double segmentMean(final int segmentIndex) {
        return get(CopyRatioParameter.SEGMENT_MEANS, CopyRatioState.SegmentMeans.class).get(segmentIndex);
    }

    boolean outlierIndicator(final int copyRatioIndex) {
        return get(CopyRatioParameter.OUTLIER_INDICATORS, CopyRatioState.OutlierIndicators.class).get(copyRatioIndex);
    }
}
