package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.model;

import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.*;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * Represents a segmented model for copy ratio fit to denoised log2 copy-ratio data.
 * The log2 copy ratios in each segment are fit by a mixture model with a normal-distribution component
 * and a uniform outlier component.  The variance of the normal-distribution component and the relative
 * contribution of the uniform outlier component in all segments are both assumed to be global parameters.
 * The mean of the normal-distribution component in each segment is taken to be a segment-level parameter.
 * The component (i.e., normal or outlier) that each copy-ratio point is drawn from is determined by a latent
 * point-level indicator.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CopyRatioModeller {
    private static final double EPSILON = 1E-10;
    private static final double LOG2_COPY_RATIO_MIN = -100.;
    private static final double LOG2_COPY_RATIO_MAX = 100.;
    private static final double LOG2_COPY_RATIO_RANGE = LOG2_COPY_RATIO_MAX - LOG2_COPY_RATIO_MIN;
    private static final double VARIANCE_MIN = EPSILON;

    private static final double OUTLIER_PROBABILITY_INITIAL = 0.05;
    private static final double OUTLIER_PROBABILITY_PRIOR_ALPHA = 5.;
    private static final double OUTLIER_PROBABILITY_PRIOR_BETA = 95.;

    private final CopyRatioData data;
    private final ParameterizedModel<CopyRatioParameter, CopyRatioState, CopyRatioData> model;

    private final List<Double> varianceSamples = new ArrayList<>();
    private final List<Double> outlierProbabilitySamples = new ArrayList<>();
    private final List<CopyRatioState.SegmentMeans> segmentMeansSamples = new ArrayList<>();
    private final List<CopyRatioState.OutlierIndicators> outlierIndicatorsSamples = new ArrayList<>();

    /**
     * Constructs a copy-ratio model given {@link CopyRatioData} with copy ratios and segments.
     * Initial point estimates of parameters are set to empirical estimates where available.
     */
    CopyRatioModeller(final CopyRatioData data) {
        this.data = Utils.nonNull(data);

        //set widths for slice sampling of variance and segment-mean posteriors using empirical variance estimate.
        //variance posterior is inverse chi-squared, segment-mean posteriors are Gaussian; the below expressions
        //approximate the standard deviations of these distributions.
        final double dataRangeOrNaN = data.getMaxLog2CopyRatioValue() - data.getMinLog2CopyRatioValue();
        final double dataRange = Double.isNaN(dataRangeOrNaN) ? LOG2_COPY_RATIO_RANGE : dataRangeOrNaN;
        final double varianceEstimateOrNaN = data.estimateVariance();
        final double varianceEstimate = Double.isNaN(varianceEstimateOrNaN) ? VARIANCE_MIN : Math.max(varianceEstimateOrNaN, VARIANCE_MIN);
        final double varianceSliceSamplingWidth = 2. * varianceEstimate;
        final double varianceMax = Math.max(10. * varianceEstimate, dataRange * dataRange);
        final double meanSliceSamplingWidth = Math.sqrt(varianceEstimate * data.getNumSegments() / data.getNumPoints());

        //the uniform log-likelihood for outliers is determined by the minimum and maximum coverages in the dataset;
        //the outlier-probability parameter should be interpreted accordingly
        final double outlierUniformLogLikelihood = -Math.log(dataRange);

        //use empirical segment means and empirical average variance across segments to initialize CopyRatioState
        final CopyRatioState initialState = new CopyRatioState(varianceEstimate, CopyRatioModeller.OUTLIER_PROBABILITY_INITIAL,
                        data.estimateSegmentMeans(), new CopyRatioState.OutlierIndicators(Collections.nCopies(data.getNumPoints(), false)));

        //define ParameterSamplers
        final ParameterSampler<Double, CopyRatioParameter, CopyRatioState, CopyRatioData> varianceSampler =
                new CopyRatioSamplers.VarianceSampler(VARIANCE_MIN, varianceMax, varianceSliceSamplingWidth);
        final ParameterSampler<Double, CopyRatioParameter, CopyRatioState, CopyRatioData> outlierProbabilitySampler =
                new CopyRatioSamplers.OutlierProbabilitySampler(OUTLIER_PROBABILITY_PRIOR_ALPHA, OUTLIER_PROBABILITY_PRIOR_BETA);
        final ParameterSampler<CopyRatioState.SegmentMeans, CopyRatioParameter, CopyRatioState, CopyRatioData> segmentMeansSampler =
                new CopyRatioSamplers.SegmentMeansSampler(LOG2_COPY_RATIO_MIN, LOG2_COPY_RATIO_MAX, meanSliceSamplingWidth);
        final ParameterSampler<CopyRatioState.OutlierIndicators, CopyRatioParameter, CopyRatioState, CopyRatioData> outlierIndicatorsSampler =
                new CopyRatioSamplers.OutlierIndicatorsSampler(outlierUniformLogLikelihood);

        model = new ParameterizedModel.GibbsBuilder<>(initialState, data)
                .addParameterSampler(CopyRatioParameter.VARIANCE, varianceSampler, Double.class)
                .addParameterSampler(CopyRatioParameter.OUTLIER_PROBABILITY, outlierProbabilitySampler, Double.class)
                .addParameterSampler(CopyRatioParameter.SEGMENT_MEANS, segmentMeansSampler, CopyRatioState.SegmentMeans.class)
                .addParameterSampler(CopyRatioParameter.OUTLIER_INDICATORS, outlierIndicatorsSampler, CopyRatioState.OutlierIndicators.class)
                .build();
    }

    /**
     * Adds {@code numSamples - numBurnIn} Markov-Chain Monte-Carlo samples of the parameter posteriors (generated using
     * Gibbs sampling) to the collections held internally.  The current {@link CopyRatioState} held internally is used
     * to initialize the Markov Chain.
     * @param numSamples    total number of samples per posterior
     * @param numBurnIn     number of burn-in samples to discard
     */
    void fitMCMC(final int numSamples,
                 final int numBurnIn) {
        ParamUtils.isPositiveOrZero(numBurnIn, "Number of burn-in samples must be non-negative.");
        Utils.validateArg(numBurnIn < numSamples, "Number of samples must be greater than number of burn-in samples.");

        //run MCMC
        final GibbsSampler<CopyRatioParameter, CopyRatioState, CopyRatioData> gibbsSampler = new GibbsSampler<>(numSamples, model);
        gibbsSampler.runMCMC();

        //update posterior samples
        varianceSamples.addAll(gibbsSampler.getSamples(CopyRatioParameter.VARIANCE, Double.class, numBurnIn));
        outlierProbabilitySamples.addAll(gibbsSampler.getSamples(CopyRatioParameter.OUTLIER_PROBABILITY, Double.class, numBurnIn));
        segmentMeansSamples.addAll(gibbsSampler.getSamples(CopyRatioParameter.SEGMENT_MEANS, CopyRatioState.SegmentMeans.class, numBurnIn));
        outlierIndicatorsSamples.addAll(gibbsSampler.getSamples(CopyRatioParameter.OUTLIER_INDICATORS, CopyRatioState.OutlierIndicators.class, numBurnIn));
    }

    /**
     * Returns an unmodifiable view of the list of samples of the variance posterior.
     * @return  unmodifiable view of the list of samples of the variance posterior
     */
    public List<Double> getVarianceSamples() {
        return Collections.unmodifiableList(varianceSamples);
    }

    /**
     * Returns an unmodifiable view of the list of samples of the outlier-probability posterior.
     * @return  unmodifiable view of the list of samples of the outlier-probability posterior
     */
    public List<Double> getOutlierProbabilitySamples() {
        return Collections.unmodifiableList(outlierProbabilitySamples);
    }

    /**
     * Returns an unmodifiable view of the list of samples of the segment-means posterior, represented as a list of
     * {@link CopyRatioState.SegmentMeans} objects.
     * @return  unmodifiable view of the list of samples of the segment-means posterior
     */
    public List<CopyRatioState.SegmentMeans> getSegmentMeansSamples() {
        return Collections.unmodifiableList(segmentMeansSamples);
    }

    /**
     * Returns an unmodifiable view of the list of samples of the outlier-indicators posterior, represented as a list of
     * {@link CopyRatioState.OutlierIndicators} objects.
     * @return  unmodifiable view of the list of samples of the outlier-indicators posterior
     */
    public List<CopyRatioState.OutlierIndicators> getOutlierIndicatorsSamples() {
        return Collections.unmodifiableList(outlierIndicatorsSamples);
    }

    /**
     * Returns a list of {@link PosteriorSummary} elements summarizing the segment-mean posterior for each segment.
     * Should only be called after {@link CopyRatioModeller#fitMCMC(int, int)} has been called.
     * @param credibleIntervalAlpha credible-interval alpha, must be in (0, 1)
     * @param ctx                   {@link JavaSparkContext} used for mllib kernel density estimation
     * @return                      list of {@link PosteriorSummary} elements summarizing the
     *                              segment-mean posterior for each segment
     */
    public List<PosteriorSummary> getSegmentMeansPosteriorSummaries(final double credibleIntervalAlpha,
                                                                    final JavaSparkContext ctx) {
        ParamUtils.inRange(credibleIntervalAlpha, 0., 1., "Credible-interval alpha must be in [0, 1].");
        Utils.nonNull(ctx);
        final int numSegments = data.getNumSegments();
        final List<PosteriorSummary> posteriorSummaries = new ArrayList<>(numSegments);
        for (int segment = 0; segment < numSegments; segment++) {
            final int j = segment;
            final List<Double> meanSamples =
                    segmentMeansSamples.stream().map(s -> s.get(j)).collect(Collectors.toList());
            posteriorSummaries.add(PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(meanSamples, credibleIntervalAlpha, ctx));
        }
        return posteriorSummaries;
    }

    /**
     * Returns a Map of {@link PosteriorSummary} elements summarizing the global parameters.
     * Should only be called after {@link CopyRatioModeller#fitMCMC(int, int)} has been called.
     * @param credibleIntervalAlpha credible-interval alpha, must be in (0, 1)
     * @param ctx                   {@link JavaSparkContext} used for mllib kernel density estimation
     * @return                      list of {@link PosteriorSummary} elements summarizing the global parameters
     */
    public Map<CopyRatioParameter, PosteriorSummary> getGlobalParameterPosteriorSummaries(final double credibleIntervalAlpha,
                                                                                          final JavaSparkContext ctx) {
        ParamUtils.inRange(credibleIntervalAlpha, 0., 1., "Credible-interval alpha must be in [0, 1].");
        Utils.nonNull(ctx);
        final Map<CopyRatioParameter, PosteriorSummary> posteriorSummaries = new LinkedHashMap<>();
        posteriorSummaries.put(CopyRatioParameter.VARIANCE, PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(varianceSamples, credibleIntervalAlpha, ctx));
        posteriorSummaries.put(CopyRatioParameter.OUTLIER_PROBABILITY, PosteriorSummaryUtils.calculateHighestPosteriorDensityAndDecilesSummary(outlierProbabilitySamples, credibleIntervalAlpha, ctx));
        return posteriorSummaries;
    }
}
