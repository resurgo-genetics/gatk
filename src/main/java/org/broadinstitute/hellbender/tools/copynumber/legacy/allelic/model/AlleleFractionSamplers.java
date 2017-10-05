package org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.model;

import org.apache.commons.math3.distribution.BetaDistribution;
import org.apache.commons.math3.random.RandomGenerator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.model.CopyRatioSegmentedData;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.mcmc.ParameterSampler;
import org.broadinstitute.hellbender.utils.mcmc.SliceSampler;

import java.util.ArrayList;
import java.util.List;
import java.util.function.Function;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

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
            return new SliceSampler(rng,
                    x -> AlleleFractionLikelihoods.logLikelihood(state.globalParameters().copyWithNewBiasVariance(x), state.minorFractions(), data),
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
            return new SliceSampler(rng,
                    x -> AlleleFractionLikelihoods.logLikelihood(state.globalParameters().copyWithNewOutlierProbability(x), state.minorFractions(), data),
                    MIN_OUTLIER_PROBABILITY, MAX_OUTLIER_PROBABILITY, outlierProbabilitySliceSamplingWidth)
                    .sample(state.outlierProbability());
        }
    }

    // sample minor fractions of all segments
    protected static final class MinorFractionsSampler implements ParameterSampler<AlleleFractionState.MinorFractions, AlleleFractionParameter, AlleleFractionState, AlleleFractionSegmentedData> {
        private static double MIN_MINOR_FRACTION = 0.;
        private static double MAX_MINOR_FRACTION = 0.5;
        private static final double PRIOR_BETA = 1.;
        private static final int NUM_POINTS_SUBSAMPLE_THRESHOLD = 1000;

        private final Function<Double, Double> logPrior;
        private final List<Double> sliceSamplingWidths;

        MinorFractionsSampler(final AlleleFractionPrior prior,
                              final List<Double> sliceSamplingWidths) {
            logPrior = f -> new BetaDistribution(null, prior.getMinorAlleleFractionPriorAlpha(), PRIOR_BETA).logDensity(2 * f);
            this.sliceSamplingWidths = sliceSamplingWidths;
        }

        @Override
        public AlleleFractionState.MinorFractions sample(final RandomGenerator rng, final AlleleFractionState state, final AlleleFractionSegmentedData data) {
            final List<Double> minorFractions = new ArrayList<>(data.getNumSegments());
            for (int segment = 0; segment < data.getNumSegments(); segment++) {
                logger.debug(String.format("Sampling minor fraction for segment %d...", segment));
                if (data.getIndexedAllelicCountsInSegment(segment).isEmpty()){
                    minorFractions.add(Double.NaN);
                } else {
                    //subsample the data if we are above the threshold
                    final List<AlleleFractionSegmentedData.IndexedAllelicCount> allelicCounts = data.getIndexedAllelicCountsInSegment(segment);
                    final List<AlleleFractionSegmentedData.IndexedAllelicCount> subsampledAllelicCounts = allelicCounts.size() > NUM_POINTS_SUBSAMPLE_THRESHOLD
                            ? IntStream.range(0, NUM_POINTS_SUBSAMPLE_THRESHOLD).boxed().map(i -> rng.nextInt(allelicCounts.size())).map(allelicCounts::get).collect(Collectors.toList())
                            : allelicCounts;
                    final double scalingFactor = (double) allelicCounts.size() / subsampledAllelicCounts.size();
                    final Function<Double, Double> segmentLogLikelihood = f -> scalingFactor * AlleleFractionLikelihoods.segmentLogLikelihood(state.globalParameters(), f, subsampledAllelicCounts);
                    final SliceSampler sampler = new SliceSampler(rng,
                            f -> logPrior.apply(f) + segmentLogLikelihood.apply(f),
                            MIN_MINOR_FRACTION, MAX_MINOR_FRACTION, sliceSamplingWidths.get(segment));
                    minorFractions.add(sampler.sample(state.segmentMinorFraction(segment)));
                }
            }
            return new AlleleFractionState.MinorFractions(minorFractions);
        }
    }
}
