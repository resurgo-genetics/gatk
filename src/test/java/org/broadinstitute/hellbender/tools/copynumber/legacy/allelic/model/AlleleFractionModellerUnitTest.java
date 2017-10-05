package org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.model;

import htsjdk.samtools.util.Log;
import org.apache.spark.api.java.JavaSparkContext;
import org.broadinstitute.hellbender.engine.spark.SparkContextFactory;
import org.broadinstitute.hellbender.utils.LoggingUtils;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;

import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * Tests the MCMC inference of the {@link AlleleFractionModeller}.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AlleleFractionModellerUnitTest extends BaseTest {
    private static final double CREDIBLE_INTERVAL_ALPHA = 0.05;

    /**
     * Test MCMC inference on simulated data.  Note that hyperparameter values used to generate the data should be recovered
     * along with outlier probability and minor fractions.
     */
    @Test
    public void testMCMC() {
        final double meanBiasSimulated = 1.2;
        final double biasVarianceSimulated = 0.04;
        testMCMC(meanBiasSimulated, biasVarianceSimulated, meanBiasSimulated, biasVarianceSimulated);
    }

    private void testMCMC(final double meanBiasSimulated, final double biasVarianceSimulated,
                          final double meanBiasExpected, final double biasVarianceExpected) {
        LoggingUtils.setLoggingLevel(Log.LogLevel.INFO);
        final JavaSparkContext ctx = SparkContextFactory.getTestSparkContext();

        final String sampleName = "test";
        final double minorAlleleFractionPriorAlpha = 1.;
        final AlleleFractionPrior prior = new AlleleFractionPrior(minorAlleleFractionPriorAlpha);

        final int numSamples = 150;
        final int numBurnIn = 50;

        final double averageHetsPerSegment = 50;
        final int numSegments = 100;
        final int averageDepth = 50;

        final double outlierProbability = 0.02;

        // note: the following tolerances could actually be made much smaller if we used more segments and/or
        // more hets -- most of the error is the sampling error of a finite simulated data set, not numerical error of MCMC
        final double minorFractionTolerance = 0.02;
        final double meanBiasTolerance = 0.02;
        final double biasVarianceTolerance = 0.01;
        final double outlierProbabilityTolerance = 0.02;
        final AlleleFractionSimulatedData simulatedData = new AlleleFractionSimulatedData(
                sampleName, averageHetsPerSegment, numSegments, averageDepth, meanBiasSimulated, biasVarianceSimulated, outlierProbability);

        final AlleleFractionModeller modeller = new AlleleFractionModeller(simulatedData.getData(), prior);
        modeller.fitMCMC(numSamples, numBurnIn);

        final List<Double> meanBiasSamples = modeller.getmeanBiasSamples();
        Assert.assertEquals(meanBiasSamples.size(), numSamples - numBurnIn);

        final List<Double> biasVarianceSamples = modeller.getBiasVarianceSamples();
        Assert.assertEquals(biasVarianceSamples.size(), numSamples - numBurnIn);

        final List<Double> outlierProbabilitySamples = modeller.getOutlierProbabilitySamples();
        Assert.assertEquals(outlierProbabilitySamples.size(), numSamples - numBurnIn);

        final List<AlleleFractionState.MinorFractions> minorFractionsSamples = modeller.getMinorFractionsSamples();
        Assert.assertEquals(minorFractionsSamples.size(), numSamples - numBurnIn);
        for (final AlleleFractionState.MinorFractions sample : minorFractionsSamples) {
            Assert.assertEquals(sample.size(), numSegments);
        }

        final List<List<Double>> minorFractionsSamplesBySegment = modeller.getMinorFractionSamplesBySegment();

        final double mcmcMeanBias = meanBiasSamples.stream().mapToDouble(x -> x).average().getAsDouble();
        final double mcmcBiasVariance = biasVarianceSamples.stream().mapToDouble(x -> x).average().getAsDouble();
        final double mcmcOutlierProbability = outlierProbabilitySamples.stream().mapToDouble(x -> x).average().getAsDouble();
        final List<Double> mcmcMinorFractions = minorFractionsSamplesBySegment
                .stream().map(list -> list.stream().mapToDouble(x -> x).average().getAsDouble())
                .collect(Collectors.toList());

        double totalSegmentError = 0.0;
        for (int segment = 0; segment < numSegments; segment++) {
            totalSegmentError += Math.abs(mcmcMinorFractions.get(segment) - simulatedData.getTrueState().segmentMinorFraction(segment));
        }

        Assert.assertEquals(mcmcMeanBias, meanBiasExpected, meanBiasTolerance);
        Assert.assertEquals(mcmcBiasVariance, biasVarianceExpected, biasVarianceTolerance);
        Assert.assertEquals(mcmcOutlierProbability, outlierProbability, outlierProbabilityTolerance);
        Assert.assertEquals(totalSegmentError / numSegments, 0.0, minorFractionTolerance);

        //test posterior summaries
        final Map<AlleleFractionParameter, PosteriorSummary> globalParameterPosteriorSummaries =
                modeller.getGlobalParameterPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);

        final PosteriorSummary meanBiasPosteriorSummary = globalParameterPosteriorSummaries.get(AlleleFractionParameter.MEAN_BIAS);
        final double meanBiasPosteriorCenter = meanBiasPosteriorSummary.getCenter();
        Assert.assertEquals(meanBiasPosteriorCenter, meanBiasExpected, meanBiasTolerance);

        final PosteriorSummary biasVariancePosteriorSummary = globalParameterPosteriorSummaries.get(AlleleFractionParameter.BIAS_VARIANCE);
        final double biasVariancePosteriorCenter = biasVariancePosteriorSummary.getCenter();
        Assert.assertEquals(biasVariancePosteriorCenter, biasVarianceExpected, biasVarianceTolerance);

        final PosteriorSummary outlierProbabilityPosteriorSummary = globalParameterPosteriorSummaries.get(AlleleFractionParameter.OUTLIER_PROBABILITY);
        final double outlierProbabilityPosteriorCenter = outlierProbabilityPosteriorSummary.getCenter();
        Assert.assertEquals(outlierProbabilityPosteriorCenter, outlierProbability, outlierProbabilityTolerance);

        final List<PosteriorSummary> minorAlleleFractionPosteriorSummaries =
                modeller.getMinorAlleleFractionsPosteriorSummaries(CREDIBLE_INTERVAL_ALPHA, ctx);
        final List<Double> minorFractionsPosteriorCenters = minorAlleleFractionPosteriorSummaries.stream().map(PosteriorSummary::getCenter).collect(Collectors.toList());
        double totalPosteriorCentersSegmentError = 0.0;
        for (int segment = 0; segment < numSegments; segment++) {
            totalPosteriorCentersSegmentError += Math.abs(minorFractionsPosteriorCenters.get(segment) - simulatedData.getTrueState().segmentMinorFraction(segment));
        }
        Assert.assertEquals(totalPosteriorCentersSegmentError / numSegments, 0.0, minorFractionTolerance);
    }
}