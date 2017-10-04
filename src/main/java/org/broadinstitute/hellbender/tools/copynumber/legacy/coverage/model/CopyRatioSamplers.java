package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.model;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.mcmc.ParameterSampler;
import org.broadinstitute.hellbender.utils.mcmc.SliceSampler;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.IntStream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class CopyRatioSamplers {
    private static final Logger logger = LogManager.getLogger(CopyRatioSamplers.class);

    private CopyRatioSamplers() {}

    //Calculates the exponent for a normal distribution; used in log-likelihood calculation below.
    private static double normalTerm(final double quantity, 
                                     final double mean, 
                                     final double variance) {
        return (quantity - mean) * (quantity - mean) / (2. * variance);
    }

    //samples log conditional posterior for the variance parameter, assuming uniform prior; this is given by
    //the product of Gaussian likelihoods for each non-outlier point t:
    //  log[product_{non-outlier t} variance^(-1/2) * exp(-(log2cr_t - mean_t)^2 / (2 * variance))] + constant
    //where mean_t is identical for all points in a segment
    static final class VarianceSampler implements ParameterSampler<Double, CopyRatioParameter, CopyRatioState, CopyRatioData> {
        private final double varianceMin;
        private final double varianceMax;
        private final double varianceSliceSamplingWidth;

        VarianceSampler(final double varianceMin, 
                        final double varianceMax, 
                        final double varianceSliceSamplingWidth) {
            this.varianceMin = varianceMin;
            this.varianceMax = varianceMax;
            this.varianceSliceSamplingWidth = varianceSliceSamplingWidth;
        }

        @Override
        public Double sample(final RandomGenerator rng, 
                             final CopyRatioState state, 
                             final CopyRatioData dataCollection) {
            logger.debug("Sampling variance...");
            final Function<Double, Double> logConditionalPDF = newVariance -> {
                final double gaussianLogNormalization = 0.5 * Math.log(newVariance);
                double ll = 0.;
                for (int segment = 0; segment < dataCollection.getNumSegments(); segment++) {
                    final List<CopyRatioData.IndexedCopyRatio> indexedCopyRatiosInSegment = dataCollection.getIndexedCopyRatiosInSegment(segment);
                    for (final CopyRatioData.IndexedCopyRatio indexedCopyRatio : indexedCopyRatiosInSegment) {
                        if (!state.outlierIndicator(indexedCopyRatio.getIndex())) {
                            ll -= normalTerm(indexedCopyRatio.getLog2CopyRatioValue(), state.segmentMean(segment), newVariance) + gaussianLogNormalization;
                        }
                    }
                }
                return ll;
            };
            return new SliceSampler(rng, logConditionalPDF, varianceMin, varianceMax, varianceSliceSamplingWidth).sample(state.variance());
        }
    }

    //samples log conditional posterior for the outlier-probability parameter, assuming Beta(alpha, beta) prior;
    //this is given by:
    //  log Beta(alpha + number of outlier points, beta + number of non-outlier points) + constant
    static final class OutlierProbabilitySampler implements ParameterSampler<Double, CopyRatioParameter, CopyRatioState, CopyRatioData> {
        private final double outlierProbabilityPriorAlpha;
        private final double outlierProbabilityPriorBeta;

        OutlierProbabilitySampler(final double outlierProbabilityPriorAlpha, 
                                  final double outlierProbabilityPriorBeta) {
            this.outlierProbabilityPriorAlpha = outlierProbabilityPriorAlpha;
            this.outlierProbabilityPriorBeta = outlierProbabilityPriorBeta;
        }

        @Override
        public Double sample(final RandomGenerator rng, 
                             final CopyRatioState state, 
                             final CopyRatioData dataCollection) {
            logger.debug("Sampling outlier probability...");
            final int numOutliers = (int) IntStream.range(0, dataCollection.getNumPoints()).filter(state::outlierIndicator).count();
            return new BetaDistribution(rng,
                    outlierProbabilityPriorAlpha + numOutliers,
                    outlierProbabilityPriorBeta + dataCollection.getNumPoints() - numOutliers).sample();
        }
    }

    //samples log conditional posteriors for the segment-mean parameters, assuming uniform priors bounded by minimum and maximum log2 copy-ratio values;
    //for each segment s, this is given by the product of Gaussian likelihoods for each non-outlier point t:
    //  log[product_{non-outlier t in s} exp(-(log2cr_t - mean_s)^2 / (2 * variance))] + constant
    static final class SegmentMeansSampler implements ParameterSampler<CopyRatioState.SegmentMeans, CopyRatioParameter, CopyRatioState, CopyRatioData> {
        private final double meanMin;
        private final double meanMax;
        private final double meanSliceSamplingWidth;

        SegmentMeansSampler(final double meanMin, 
                            final double meanMax, 
                            final double meanSliceSamplingWidth) {
            this.meanMin = meanMin;
            this.meanMax = meanMax;
            this.meanSliceSamplingWidth = meanSliceSamplingWidth;
        }

        @Override
        public CopyRatioState.SegmentMeans sample(final RandomGenerator rng,
                                                  final CopyRatioState state, 
                                                  final CopyRatioData dataCollection) {
            final List<Double> means = new ArrayList<>(dataCollection.getNumSegments());
            for (int segment = 0; segment < dataCollection.getNumSegments(); segment++) {
                final List<CopyRatioData.IndexedCopyRatio> indexedCopyRatiosInSegment = dataCollection.getIndexedCopyRatiosInSegment(segment);
                if (indexedCopyRatiosInSegment.isEmpty()) {
                    means.add(Double.NaN);
                } else {
                    logger.debug(String.format("Sampling mean for segment %d...", segment));
                    final Function<Double, Double> logConditionalPDF = newMean ->
                            indexedCopyRatiosInSegment.stream()
                                    .filter(c -> !state.outlierIndicator(c.getIndex()))
                                    .mapToDouble(c -> -normalTerm(c.getLog2CopyRatioValue(), newMean, state.variance()))
                                    .sum();
                    final SliceSampler sampler = new SliceSampler(rng, logConditionalPDF, meanMin, meanMax, meanSliceSamplingWidth);
                    means.add(sampler.sample(state.segmentMean(segment)));
                }
            }
            return new CopyRatioState.SegmentMeans(means);
        }
    }

    //samples log conditional posteriors for the outlier-indicator parameters; for each point t, this is given by:
    //          z_t * [log outlier_prob + outlierUniformLogLikelihood]
    //  + (1 - z_t) * [log(1 - outlier_prob) - log(2 * pi * variance)/2 - (log2cr_t - mean_t)^2 / (2 * variance)]
    //  + const
    //where z_t is the indicator for point t, and outlier_prob is the outlier probability.
    //note that we compute the normalizing constant, so that we can sample a new indicator value by simply sampling
    //uniformly in [0, 1] and checking whether the resulting value is less than the probability of being an outlier
    //(corresponding to the first line in the unnormalized expression above)
    static final class OutlierIndicatorsSampler implements ParameterSampler<CopyRatioState.OutlierIndicators, CopyRatioParameter, CopyRatioState, CopyRatioData> {
        private final double outlierUniformLogLikelihood;

        OutlierIndicatorsSampler(final double outlierUniformLogLikelihood) {
            this.outlierUniformLogLikelihood = outlierUniformLogLikelihood;
        }

        @Override
        public CopyRatioState.OutlierIndicators sample(final RandomGenerator rng,
                                                       final CopyRatioState state,
                                                       final CopyRatioData dataCollection) {
            logger.debug("Sampling outlier indicators...");
            final double outlierUnnormalizedLogProbability =
                    Math.log(state.outlierProbability()) + outlierUniformLogLikelihood;
            final double notOutlierUnnormalizedLogProbabilityPrefactor =
                    Math.log(1. - state.outlierProbability()) - 0.5 * Math.log(2 * Math.PI * state.variance());
            final List<Boolean> indicators = new ArrayList<>();
            for (int segment = 0; segment < dataCollection.getNumSegments(); segment++) {
                final List<CopyRatioData.IndexedCopyRatio> indexedCopyRatiosInSegment = dataCollection.getIndexedCopyRatiosInSegment(segment);
                for (final CopyRatioData.IndexedCopyRatio indexedCopyRatio : indexedCopyRatiosInSegment) {
                    final double notOutlierUnnormalizedLogProbability =
                            notOutlierUnnormalizedLogProbabilityPrefactor
                                    - normalTerm(indexedCopyRatio.getLog2CopyRatioValue(), state.segmentMean(segment), state.variance());
                    //note: we are working in natural log space, so we divide by ln(10) before using normalizeFromLog10
                    final double conditionalProbability =
                            MathUtils.normalizeFromLog10ToLinearSpace(new double[]{
                                    MathUtils.logToLog10(outlierUnnormalizedLogProbability),
                                    MathUtils.logToLog10(notOutlierUnnormalizedLogProbability)})[0];
                    indicators.add(rng.nextDouble() < conditionalProbability);
                }
            }
            return new CopyRatioState.OutlierIndicators(indicators);
        }
    }
}