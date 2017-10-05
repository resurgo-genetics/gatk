package org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.model;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.mcmc.ParameterSampler;
import org.broadinstitute.hellbender.utils.mcmc.SliceSampler;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Sampler classes for the allele-fraction model.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class AlleleFractionSamplers {
    private static final Logger logger = LogManager.getLogger(AlleleFractionSamplers.class);

    private AlleleFractionSamplers() {}

    static final class MeanBiasSampler implements ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionSegmentedData> {
        private static final double MIN_MEAN_BIAS = 0.;

        private final double maxMeanBias;
        private final double meanBiasSliceSamplingWidth;

        MeanBiasSampler(final double maxMeanBias,
                        final double meanBiasSliceSamplingWidth) {
            this.maxMeanBias = maxMeanBias;
            this.meanBiasSliceSamplingWidth = meanBiasSliceSamplingWidth;
        }

        @Override
        public Double sample(final RandomGenerator rng,
                             final AlleleFractionState state,
                             final AlleleFractionSegmentedData data) {
            logger.debug("Sampling mean bias...");
            return new SliceSampler(rng,
                    x -> AlleleFractionLikelihoods.logLikelihood(state.globalParameters().copyWithNewMeanBias(x),  state.minorFractions(), data),
                    MIN_MEAN_BIAS, maxMeanBias, meanBiasSliceSamplingWidth)
                    .sample(state.meanBias());
        }
    }

    static final class BiasVarianceSampler implements ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionSegmentedData> {
        private static final double MIN_BIAS_VARIANCE = 1E-10;

        private final double maxBiasVariance;
        private final double biasVarianceSliceSamplingWidth;

        BiasVarianceSampler(final double maxBiasVariance,
                            final double biasVarianceSliceSamplingWidth) {
            this.maxBiasVariance = maxBiasVariance;
            this.biasVarianceSliceSamplingWidth = biasVarianceSliceSamplingWidth;
        }

        @Override
        public Double sample(final RandomGenerator rng,
                             final AlleleFractionState state,
                             final AlleleFractionSegmentedData data) {
            logger.debug("Sampling bias variance...");
            return new SliceSampler(rng, x -> AlleleFractionLikelihoods.logLikelihood(
                    state.globalParameters().copyWithNewBiasVariance(x), state.minorFractions(), data),
                    MIN_BIAS_VARIANCE, maxBiasVariance, biasVarianceSliceSamplingWidth)
                    .sample(state.biasVariance());
        }
    }

    protected static final class OutlierProbabilitySampler implements ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionSegmentedData> {
        private static final double MIN_OUTLIER_PROBABILITY = 0.;
        private static final double MAX_OUTLIER_PROBABILITY = 1.;

        private final double outlierProbabilitySliceSamplingWidth;

        OutlierProbabilitySampler(final double outlierProbabilitySliceSamplingWidth) {
            this.outlierProbabilitySliceSamplingWidth = outlierProbabilitySliceSamplingWidth;
        }

        @Override
        public Double sample(final RandomGenerator rng,
                             final AlleleFractionState state,
                             final AlleleFractionSegmentedData data) {
            logger.debug("Sampling outlier probability...");
            return new SliceSampler(rng, x -> AlleleFractionLikelihoods.logLikelihood(
                    state.globalParameters().copyWithNewOutlierProbability(x), state.minorFractions(), data),
                    MIN_OUTLIER_PROBABILITY, MAX_OUTLIER_PROBABILITY, outlierProbabilitySliceSamplingWidth)
                    .sample(state.outlierProbability());
        }
    }

    // sample minor fractions of all segments
    protected static final class MinorFractionsSampler implements ParameterSampler<AlleleFractionState.MinorFractions, AlleleFractionParameter, AlleleFractionState, AlleleFractionSegmentedData> {
        private final List<PerSegmentMinorFractionSampler> perSegmentSamplers = new ArrayList<>();

        MinorFractionsSampler(final AlleleFractionPrior prior,
                              final List<Double> sliceSamplingWidths) {
            final int numSegments = sliceSamplingWidths.size();
            for (int segment = 0; segment < numSegments; segment++) {
                perSegmentSamplers.add(new PerSegmentMinorFractionSampler(segment, prior.getMinorAlleleFractionPriorAlpha(), sliceSamplingWidths.get(segment)));
            }
        }

        @Override
        public AlleleFractionState.MinorFractions sample(final RandomGenerator rng, final AlleleFractionState state, final AlleleFractionSegmentedData data) {
            return new AlleleFractionState.MinorFractions(perSegmentSamplers.stream()
                    .map(sampler -> sampler.sample(rng, state, data)).collect(Collectors.toList()));
        }
    }

    // sample minor fraction of a single segment
    private static final class PerSegmentMinorFractionSampler implements ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionSegmentedData> {
        private static double MIN_MINOR_FRACTION = 0.;
        private static double MAX_MINOR_FRACTION = 0.5;
        private static final double PRIOR_BETA = 1.;

        private final int segmentIndex;
        private final double sliceSamplingWidth;
        private final Function<Double, Double> logPrior;

        PerSegmentMinorFractionSampler(final int segmentIndex,
                                       final double priorAlpha,
                                       final double sliceSamplingWidth) {
            this.segmentIndex = segmentIndex;
            this.sliceSamplingWidth = sliceSamplingWidth;
            logPrior = f -> new BetaDistribution(null, priorAlpha, PRIOR_BETA).logDensity(2 * f);
        }

        @Override
        public Double sample(final RandomGenerator rng, final AlleleFractionState state, final AlleleFractionSegmentedData data) {
            logger.debug(String.format("Sampling minor fraction for segment %d...", segmentIndex));
            if (data.getIndexedAllelicCountsInSegment(segmentIndex).isEmpty()) {
                return Double.NaN;
            }
            return new SliceSampler(rng, f -> logPrior.apply(f) + AlleleFractionLikelihoods.segmentLogLikelihood(
                    state.globalParameters(), f, segmentIndex, data),
                    MIN_MINOR_FRACTION, MAX_MINOR_FRACTION, sliceSamplingWidth)
                    .sample(state.segmentMinorFraction(segmentIndex));
        }
    }
}
