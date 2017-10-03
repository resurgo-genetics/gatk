package org.broadinstitute.hellbender.tools.exome.allelefraction;

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

    // sample mean bias
    protected static final class MeanBiasSampler implements ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionData> {
        private static final double MIN_MEAN_BIAS = 0.;

        private final double maxMeanBias;
        private final double meanBiasSliceSamplingWidth;
        private final int numPointsSubsamplingLimit;
        private AlleleFractionData dataGlobal;

        MeanBiasSampler(final double maxMeanBias, final double meanBiasSliceSamplingWidth, final int numPointsSubsamplingLimit) {
            this.maxMeanBias = maxMeanBias;
            this.meanBiasSliceSamplingWidth = meanBiasSliceSamplingWidth;
            this.numPointsSubsamplingLimit = numPointsSubsamplingLimit;
            dataGlobal = null;
        }

        @Override
        public Double sample(final RandomGenerator rng, final AlleleFractionState state, final AlleleFractionData data) {
            logger.debug("Sampling mean bias...");
            final AllelicPanelOfNormals allelicPoN = data.getPoN();
            if (allelicPoN.equals(AllelicPanelOfNormals.EMPTY_PON)) {
                dataGlobal = data.getNumPoints() >= numPointsSubsamplingLimit && dataGlobal == null
                        ? data.subsample(rng, numPointsSubsamplingLimit)
                        : data;
                final double scalingFactor = (double) numPointsSubsamplingLimit / data.getNumPoints();
                return new SliceSampler(rng, x -> scalingFactor * AlleleFractionLikelihoods.logLikelihood(
                        state.globalParameters().copyWithNewMeanBias(x),  state.minorFractions(), dataGlobal),
                        MIN_MEAN_BIAS, maxMeanBias, meanBiasSliceSamplingWidth)
                        .sample(state.meanBias());
            }
            return allelicPoN.getGlobalMeanBias(); // if PoN is available, always return MLE mean bias as "sample"
        }
    }

    // sample bias variance
    protected static final class BiasVarianceSampler implements ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionData> {
        private static final double MIN_BIAS_VARIANCE = 1E-10;

        private final double maxBiasVariance;
        private final double biasVarianceSliceSamplingWidth;
        private final int numPointsSubsamplingLimit;
        private AlleleFractionData dataGlobal;

        BiasVarianceSampler(final double maxBiasVariance, final double biasVarianceSliceSamplingWidth, final int numPointsSubsamplingLimit) {
            this.maxBiasVariance = maxBiasVariance;
            this.biasVarianceSliceSamplingWidth = biasVarianceSliceSamplingWidth;
            this.numPointsSubsamplingLimit = numPointsSubsamplingLimit;
            dataGlobal = null;
        }

        @Override
        public Double sample(final RandomGenerator rng, final AlleleFractionState state, final AlleleFractionData data) {
            logger.debug("Sampling bias variance...");
            final AllelicPanelOfNormals allelicPoN = data.getPoN();
            if (allelicPoN.equals(AllelicPanelOfNormals.EMPTY_PON)) {
                dataGlobal = data.getNumPoints() >= numPointsSubsamplingLimit && dataGlobal == null
                        ? data.subsample(rng, numPointsSubsamplingLimit)
                        : data;
                final double scalingFactor = (double) numPointsSubsamplingLimit / data.getNumPoints();
                return new SliceSampler(rng, x -> scalingFactor * AlleleFractionLikelihoods.logLikelihood(
                        state.globalParameters().copyWithNewBiasVariance(x), state.minorFractions(), dataGlobal),
                        MIN_BIAS_VARIANCE, maxBiasVariance, biasVarianceSliceSamplingWidth)
                        .sample(state.biasVariance());
            }
            return allelicPoN.getGlobalBiasVariance(); // if PoN is available, always return MLE bias variance as "sample"
        }
    }

    // sample outlier probability
    protected static final class OutlierProbabilitySampler implements ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionData> {
        private static final double MIN_OUTLIER_PROBABILITY = 0.;
        private static final double MAX_OUTLIER_PROBABILITY = 1.;

        private final double outlierProbabilitySliceSamplingWidth;
        private final int numPointsSubsamplingLimit;
        private AlleleFractionData dataGlobal;

        OutlierProbabilitySampler(final double outlierProbabilitySliceSamplingWidth, final int numPointsSubsamplingLimit) {
            this.outlierProbabilitySliceSamplingWidth = outlierProbabilitySliceSamplingWidth;
            this.numPointsSubsamplingLimit = numPointsSubsamplingLimit;
            dataGlobal = null;
        }

        @Override
        public Double sample(final RandomGenerator rng, final AlleleFractionState state, final AlleleFractionData data) {
            logger.debug("Sampling outlier probability...");
            dataGlobal = data.getNumPoints() >= numPointsSubsamplingLimit && dataGlobal == null
                    ? data.subsample(rng, numPointsSubsamplingLimit)
                    : data;
            final double scalingFactor = (double) numPointsSubsamplingLimit / data.getNumPoints();
            return new SliceSampler(rng, x -> scalingFactor * AlleleFractionLikelihoods.logLikelihood(
                    state.globalParameters().copyWithNewOutlierProbability(x), state.minorFractions(), dataGlobal),
                    MIN_OUTLIER_PROBABILITY, MAX_OUTLIER_PROBABILITY, outlierProbabilitySliceSamplingWidth)
                    .sample(state.outlierProbability());
        }
    }

    // sample minor fraction of a single segment
    private static final class PerSegmentMinorFractionSampler implements ParameterSampler<Double, AlleleFractionParameter, AlleleFractionState, AlleleFractionData> {
        private static double MIN_MINOR_FRACTION = 0.;
        private static double MAX_MINOR_FRACTION = 0.5;
        private static final double PRIOR_BETA = 1.;

        private final int segmentIndex;
        private final double sliceSamplingWidth;
        private final Function<Double, Double> logPrior;

        PerSegmentMinorFractionSampler(final int segmentIndex, final double priorAlpha, final double sliceSamplingWidth) {
            this.segmentIndex = segmentIndex;
            this.sliceSamplingWidth = sliceSamplingWidth;
            logPrior = f -> new BetaDistribution(null, priorAlpha, PRIOR_BETA).logDensity(2 * f);
        }

        @Override
        public Double sample(final RandomGenerator rng, final AlleleFractionState state, final AlleleFractionData data) {
            logger.debug(String.format("Sampling minor fraction for segment %d...", segmentIndex));
            if (data.getNumHetsInSegment(segmentIndex) == 0) {
                return Double.NaN;
            }
            return new SliceSampler(rng, f -> logPrior.apply(f) + AlleleFractionLikelihoods.segmentLogLikelihood(
                    state.globalParameters(), f, data.getCountsInSegment(segmentIndex), data.getPoN()),
                    MIN_MINOR_FRACTION, MAX_MINOR_FRACTION, sliceSamplingWidth)
                    .sample(state.segmentMinorFraction(segmentIndex));
        }
    }

    // sample minor fractions of all segments
    protected static final class MinorFractionsSampler implements ParameterSampler<AlleleFractionState.MinorFractions, AlleleFractionParameter, AlleleFractionState, AlleleFractionData> {
        private final List<PerSegmentMinorFractionSampler> perSegmentSamplers = new ArrayList<>();

        MinorFractionsSampler(final double minorAlleleFractionPriorAlpha, final List<Double> sliceSamplingWidths) {
            final int numSegments = sliceSamplingWidths.size();
            for (int segment = 0; segment < numSegments; segment++) {
                perSegmentSamplers.add(new PerSegmentMinorFractionSampler(segment, minorAlleleFractionPriorAlpha, sliceSamplingWidths.get(segment)));
            }
        }

        @Override
        public AlleleFractionState.MinorFractions sample(final RandomGenerator rng, final AlleleFractionState state, final AlleleFractionData data) {
            return new AlleleFractionState.MinorFractions(perSegmentSamplers.stream()
                    .map(sampler -> sampler.sample(rng, state, data)).collect(Collectors.toList()));
        }
    }
}
