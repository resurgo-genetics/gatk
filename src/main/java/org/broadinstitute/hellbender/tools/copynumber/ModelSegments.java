package org.broadinstitute.hellbender.tools.copynumber;

import com.google.common.collect.ImmutableSet;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.BetaFeature;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.allelic.model.AlleleFractionPrior;
import org.broadinstitute.hellbender.tools.copynumber.allelic.segmentation.AlleleFractionKernelSegmenter;
import org.broadinstitute.hellbender.tools.copynumber.allelic.segmentation.AlleleFractionSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.coverage.copyratio.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.coverage.segmentation.CopyRatioKernelSegmenter;
import org.broadinstitute.hellbender.tools.copynumber.coverage.segmentation.CopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.LocatableCollection;
import org.broadinstitute.hellbender.tools.copynumber.multidimensional.model.CRAFModeller;
import org.broadinstitute.hellbender.tools.copynumber.multidimensional.model.ModeledSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.multidimensional.segmentation.CRAFSegment;
import org.broadinstitute.hellbender.tools.copynumber.multidimensional.segmentation.CRAFSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.multidimensional.segmentation.MultidimensionalKernelSegmenter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
@CommandLineProgramProperties(
        summary = "Model segmented copy ratio from denoised read counts and minor-allele fraction from allelic counts.",
        oneLineSummary = "Model segmented copy ratio from denoised read counts.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
@BetaFeature
public final class ModelSegments extends CommandLineProgram {
    private static final double EPSILON = 1E-10;

    //filename tags for output
    public static final String HET_ALLELIC_COUNTS_FILE_SUFFIX = ".hets.tsv";
    public static final String SEGMENTS_FILE_SUFFIX = ".seg";
    public static final String COPY_RATIO_SEGMENTS_FILE_SUFFIX = ".cr" + SEGMENTS_FILE_SUFFIX;
    public static final String ALLELE_FRACTION_SEGMENTS_FILE_SUFFIX = ".af" + SEGMENTS_FILE_SUFFIX;
    public static final String CRAF_SEGMENTS_FILE_SUFFIX = ".craf" + SEGMENTS_FILE_SUFFIX;
    public static final String BEGIN_FIT_FILE_TAG = ".modelBegin";
    public static final String FINAL_FIT_FILE_TAG = ".modelFinal";
    public static final String CR_PARAMETER_FILE_SUFFIX = ".cr.param";
    public static final String AF_PARAMETER_FILE_SUFFIX = ".af.param";

    public static final String MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_LONG_NAME = "maxNumSegmentsPerChromosome";
    public static final String MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_SHORT_NAME = "maxNumSegsPerChr";

    public static final String MINIMUM_TOTAL_ALLELE_COUNT_LONG_NAME = "minTotalAlleleCount";
    public static final String MINIMUM_TOTAL_ALLELE_COUNT_SHORT_NAME = "minAC";

    public static final String GENOTYPING_P_VALUE_THRESHOLD_LONG_NAME = "genotypingPValueThreshold";
    public static final String GENOTYPING_P_VALUE_THRESHOLD_SHORT_NAME = "pValTh";

    public static final String GENOTYPING_BASE_ERROR_RATE_LONG_NAME = "genotypingBaseErrorRate";
    public static final String GENOTYPING_BASE_ERROR_RATE_SHORT_NAME = "baseErrRate";

    public static final String KERNEL_VARIANCE_COPY_RATIO_LONG_NAME = "kernelVarianceCopyRatio";
    public static final String KERNEL_VARIANCE_COPY_RATIO_SHORT_NAME = "kernVarCR";

    public static final String KERNEL_VARIANCE_ALLELE_FRACTION_LONG_NAME = "kernelVarianceAlleleFraction";
    public static final String KERNEL_VARIANCE_ALLELE_FRACTION_SHORT_NAME = "kernVarAF";

    public static final String KERNEL_APPROXIMATION_DIMENSION_LONG_NAME = "kernelApproximationDimension";
    public static final String KERNEL_APPROXIMATION_DIMENSION_SHORT_NAME = "kernApproxDim";

    public static final String WINDOW_SIZES_LONG_NAME = "windowSizes";
    public static final String WINDOW_SIZES_SHORT_NAME = "winSizes";

    public static final String NUM_CHANGEPOINTS_PENALTY_FACTOR_COPY_RATIO_LONG_NAME = "numChangepointsPenaltyFactorCopyRatio";
    public static final String NUM_CHANGEPOINTS_PENALTY_FACTOR_COPY_RATIO_SHORT_NAME = "numChangeptsPenCR";

    public static final String NUM_CHANGEPOINTS_PENALTY_FACTOR_ALLELE_FRACTION_LONG_NAME = "numChangepointsPenaltyFactorAlleleFraction";
    public static final String NUM_CHANGEPOINTS_PENALTY_FACTOR_ALLELE_FRACTION_SHORT_NAME = "numChangeptsPenAF";

    public static final String MINOR_ALLELE_FRACTION_PRIOR_ALPHA_LONG_NAME = "minorAlleleFractionPriorAlpha";
    public static final String MINOR_ALLELE_FRACTION_PRIOR_ALPHA_SHORT_NAME = "alphaAF";

    public static final String NUM_SAMPLES_COPY_RATIO_LONG_NAME = "numSamplesCopyRatio";
    public static final String NUM_SAMPLES_COPY_RATIO_SHORT_NAME = "numSampCR";

    public static final String NUM_BURN_IN_COPY_RATIO_LONG_NAME = "numBurnInCopyRatio";
    public static final String NUM_BURN_IN_COPY_RATIO_SHORT_NAME = "numBurnCR";

    public static final String NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME = "numSamplesAlleleFraction";
    public static final String NUM_SAMPLES_ALLELE_FRACTION_SHORT_NAME = "numSampAF";

    public static final String NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME = "numBurnInAlleleFraction";
    public static final String NUM_BURN_IN_ALLELE_FRACTION_SHORT_NAME = "numBurnAF";

    public static final String SMOOTHING_CREDIBLE_INTERVAL_THRESHOLD_COPY_RATIO_LONG_NAME = "smoothingThresholdCopyRatio";
    public static final String SMOOTHING_CREDIBLE_INTERVAL_THRESHOLD_COPY_RATIO_SHORT_NAME = "smoothThCR";

    public static final String SMOOTHING_CREDIBLE_INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME = "smoothingThresholdAlleleFraction";
    public static final String SMOOTHING_CREDIBLE_INTERVAL_THRESHOLD_ALLELE_FRACTION_SHORT_NAME = "smoothThAF";

    public static final String MAX_NUM_SMOOTHING_ITERATIONS_LONG_NAME = "maxNumSmoothingIterations";
    public static final String MAX_NUM_SMOOTHING_ITERATIONS_SHORT_NAME = "maxNumSmoothIter";

    public static final String NUM_SMOOTHING_ITERATIONS_PER_FIT_LONG_NAME = "numSmoothingIterationsPerFit";
    public static final String NUM_SMOOTHING_ITERATIONS_PER_FIT_SHORT_NAME = "numSmoothIterPerFit";

    @Argument(
            doc = "Input file containing denoised copy-ratio profile (output of DenoiseReadCounts).",
            fullName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_LONG_NAME,
            shortName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME,
            optional = true
    )
    private File inputDenoisedCopyRatiosFile = null;

    @Argument(
            doc = "Input file containing allelic counts (output of CollectAllelicCounts).",
            fullName = CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = true
    )
    private File inputAllelicCountsFile = null;

    @Argument(
            doc = "Input file containing allelic counts for a matched normal (output of CollectAllelicCounts).",
            fullName = CopyNumberStandardArgument.NORMAL_ALLELIC_COUNTS_FILE_LONG_NAME,
            shortName = CopyNumberStandardArgument.NORMAL_ALLELIC_COUNTS_FILE_SHORT_NAME,
            optional = true
    )
    private File inputNormalAllelicCountsFile = null;

    @Argument(
            doc = "Prefix for output files.",
            fullName =  CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME,
            shortName = CopyNumberStandardArgument.OUTPUT_PREFIX_SHORT_NAME
    )
    private String outputPrefix;

    @Argument(
            doc = "Output directory.",
            fullName =  StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputDir;

    @Argument(
            doc = "Maximum number of segments allowed per chromosome.",
            fullName = MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_LONG_NAME,
            shortName = MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_SHORT_NAME,
            minValue = 1,
            optional = true
    )
    private int maxNumSegmentsPerChromosome = 1000;

    @Argument(
            doc = "Minimum total count for filtering allelic counts.",
            fullName = MINIMUM_TOTAL_ALLELE_COUNT_LONG_NAME,
            shortName = MINIMUM_TOTAL_ALLELE_COUNT_SHORT_NAME,
            minValue = 0,
            optional = true
    )
    private int minTotalAlleleCount = 30;

    @Argument(
            doc = "P-value threshold for genotyping and filtering homozygous allelic counts.",
            fullName = GENOTYPING_P_VALUE_THRESHOLD_LONG_NAME,
            shortName = GENOTYPING_P_VALUE_THRESHOLD_SHORT_NAME,
            optional = true
    )
    private double genotypingPValueThreshold = 1E-2;

    @Argument(
            doc = "Base error rate for genotyping and filtering homozygous allelic counts.",
            fullName = GENOTYPING_BASE_ERROR_RATE_LONG_NAME,
            shortName = GENOTYPING_BASE_ERROR_RATE_SHORT_NAME,
            optional = true
    )
    private double genotypingBaseErrorRate = 1E-2;

    @Argument(
            doc = "Variance of Gaussian kernel for copy-ratio segmentation.  If zero, a linear kernel will be used.",
            fullName = KERNEL_VARIANCE_COPY_RATIO_LONG_NAME,
            shortName = KERNEL_VARIANCE_COPY_RATIO_SHORT_NAME,
            minValue = 0.,
            optional = true
    )
    private double kernelVarianceCopyRatio = 0.;

    @Argument(
            doc = "Variance of Gaussian kernel for allele-fraction segmentation.  If zero, a linear kernel will be used.",
            fullName = KERNEL_VARIANCE_ALLELE_FRACTION_LONG_NAME,
            shortName = KERNEL_VARIANCE_ALLELE_FRACTION_SHORT_NAME,
            minValue = 0.,
            optional = true
    )
    private double kernelVarianceAlleleFraction = 0.01;

    @Argument(
            doc = "Dimension of kernel approximation.  A subsample containing this number of data points " +
                    "will be used to construct the approximation for each chromosome.  " +
                    "If the total number of data points in a chromosome is greater " +
                    "than this number, then all data points in the chromosome will be used.  " +
                    "Time complexity scales quadratically and space complexity scales linearly with this parameter.",
            fullName = KERNEL_APPROXIMATION_DIMENSION_LONG_NAME,
            shortName = KERNEL_APPROXIMATION_DIMENSION_SHORT_NAME,
            minValue = 1,
            optional = true
    )
    private int kernelApproximationDimension = 100;

    @Argument(
            doc = "Window sizes to use for calculating local changepoint costs.  " +
                    "For each window size, the cost for each data point to be a changepoint will be calculated " +
                    "assuming that it demarcates two adjacent segments of that size.  " +
                    "Including small (large) window sizes will increase sensitivity to small (large) events.  " +
                    "Duplicate values will be ignored.",
            fullName = WINDOW_SIZES_LONG_NAME,
            shortName = WINDOW_SIZES_SHORT_NAME,
            minValue = 1,
            optional = true
    )
    private List<Integer> windowSizes = new ArrayList<>(Arrays.asList(8, 16, 32, 64, 128, 256));

    @Argument(
            doc = "Factor A for the penalty on the number of changepoints per chromosome for copy-ratio segmentation.  " +
                    "Adds a penalty of the form A * C * [1 + log (N / C)], " +
                    "where C is the number of changepoints in the chromosome, " +
                    "to the cost function for each chromosome.  " +
                    "Must be non-negative.",
            fullName = NUM_CHANGEPOINTS_PENALTY_FACTOR_COPY_RATIO_LONG_NAME,
            shortName = NUM_CHANGEPOINTS_PENALTY_FACTOR_COPY_RATIO_SHORT_NAME,
            minValue = 0.,
            optional = true
    )
    private double numChangepointsPenaltyFactorCopyRatio = 1.;

    @Argument(
            doc = "Factor A for the penalty on the number of changepoints per chromosome for allele-fraction segmentation.  " +
                    "Adds a penalty of the form A * C * [1 + log (N / C)], " +
                    "where C is the number of changepoints in the chromosome, " +
                    "to the cost function for each chromosome.  " +
                    "Must be non-negative.",
            fullName = NUM_CHANGEPOINTS_PENALTY_FACTOR_ALLELE_FRACTION_LONG_NAME,
            shortName = NUM_CHANGEPOINTS_PENALTY_FACTOR_ALLELE_FRACTION_SHORT_NAME,
            minValue = 0.,
            optional = true
    )
    private double numChangepointsPenaltyFactorAlleleFraction = 1.;

    @Argument(
            doc = "Alpha hyperparameter for the 4-parameter beta-distribution prior on segment minor-allele fraction. " +
                    "The prior for the minor-allele fraction f in each segment is assumed to be Beta(alpha, 1, 0, 1/2). " +
                    "Increasing this hyperparameter will reduce the effect of reference bias at the expense of sensitivity.",
            fullName = MINOR_ALLELE_FRACTION_PRIOR_ALPHA_LONG_NAME,
            shortName = MINOR_ALLELE_FRACTION_PRIOR_ALPHA_SHORT_NAME,
            optional = true,
            minValue = 1
    )
    private double minorAlleleFractionPriorAlpha = 25.;

    @Argument(
            doc = "Total number of MCMC samples for copy-ratio model.",
            fullName = NUM_SAMPLES_COPY_RATIO_LONG_NAME,
            shortName = NUM_SAMPLES_COPY_RATIO_SHORT_NAME,
            optional = true,
            minValue = 1
    )
    private int numSamplesCopyRatio = 100;

    @Argument(
            doc = "Number of burn-in samples to discard for copy-ratio model.",
            fullName = NUM_BURN_IN_COPY_RATIO_LONG_NAME,
            shortName = NUM_BURN_IN_COPY_RATIO_SHORT_NAME,
            optional = true,
            minValue = 0
    )
    private int numBurnInCopyRatio = 50;

    @Argument(
            doc = "Total number of MCMC samples for allele-fraction model.",
            fullName = NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME,
            shortName = NUM_SAMPLES_ALLELE_FRACTION_SHORT_NAME,
            optional = true,
            minValue = 1
    )
    private int numSamplesAlleleFraction = 100;

    @Argument(
            doc = "Number of burn-in samples to discard for allele-fraction model.",
            fullName = NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME,
            shortName = NUM_BURN_IN_ALLELE_FRACTION_SHORT_NAME,
            optional = true,
            minValue = 0
    )
    private int numBurnInAlleleFraction = 50;

    @Argument(
            doc = "Number of 10% equal-tailed credible-interval widths to use for copy-ratio segmentation smoothing.",
            fullName = SMOOTHING_CREDIBLE_INTERVAL_THRESHOLD_COPY_RATIO_LONG_NAME,
            shortName = SMOOTHING_CREDIBLE_INTERVAL_THRESHOLD_COPY_RATIO_SHORT_NAME,
            optional = true,
            minValue = 0.
    )
    private double smoothingCredibleIntervalThresholdCopyRatio = 4.;

    @Argument(
            doc = "Number of 10% equal-tailed credible-interval widths to use for allele-fraction segmentation smoothing.",
            fullName = SMOOTHING_CREDIBLE_INTERVAL_THRESHOLD_ALLELE_FRACTION_LONG_NAME,
            shortName = SMOOTHING_CREDIBLE_INTERVAL_THRESHOLD_ALLELE_FRACTION_SHORT_NAME,
            optional = true,
            minValue = 0.
    )
    private double smoothingCredibleIntervalThresholdAlleleFraction = 2.;

    @Argument(
            doc = "Maximum number of iterations allowed for segmentation smoothing.",
            fullName = MAX_NUM_SMOOTHING_ITERATIONS_LONG_NAME,
            shortName = MAX_NUM_SMOOTHING_ITERATIONS_SHORT_NAME,
            optional = true,
            minValue = 0
    )
    private int maxNumSmoothingIterations = 10;

    @Argument(
            doc = "Number of segmentation-smoothing iterations per MCMC model refit. " +
                    "(Increasing this will decrease runtime, but the final number of segments may be higher. " +
                    "Setting this to 0 will completely disable model refitting between iterations.)",
            fullName = NUM_SMOOTHING_ITERATIONS_PER_FIT_LONG_NAME,
            shortName = NUM_SMOOTHING_ITERATIONS_PER_FIT_SHORT_NAME,
            optional = true,
            minValue = 0
    )
    private int numSmoothingIterationsPerFit = 0;

    //initialize data/segment variables, some of which may be optional
    private CopyRatioCollection denoisedCopyRatios = null;
    private AllelicCountCollection hetAllelicCounts = null;

    @Override
    protected Object doWork() {
        validateArguments();

        //perform one-dimensional or multidimensional segmentation as appropriate and write the combined segments to file
        final CRAFSegmentCollection crafSegments;
        if (inputDenoisedCopyRatiosFile != null && inputAllelicCountsFile == null) {
            readDenoisedCopyRatios();
            final CopyRatioSegmentCollection copyRatioSegments = performCopyRatioSegmentation();
            crafSegments = new CRAFSegmentCollection(
                    copyRatioSegments.getSampleName(),
                    copyRatioSegments.getRecords().stream()
                            .map(s -> new CRAFSegment(s.getInterval(), s.getNumPoints(), 0, s.getMeanLog2CopyRatio()))
                            .collect(Collectors.toList()));
            hetAllelicCounts = new AllelicCountCollection(denoisedCopyRatios.getSampleName(), Collections.emptyList()); //create an empty collection with the appropriate name
        } else if (inputDenoisedCopyRatiosFile == null && inputAllelicCountsFile != null) {
            readAndFilterAllelicCounts();
            final AlleleFractionSegmentCollection alleleFractionSegments = performAlleleFractionSegmentation();
            crafSegments = new CRAFSegmentCollection(
                    alleleFractionSegments.getSampleName(),
                    alleleFractionSegments.getRecords().stream()
                            .map(s -> new CRAFSegment(s.getInterval(), 0, s.getNumPoints(), Double.NaN))
                            .collect(Collectors.toList()));
            denoisedCopyRatios = new CopyRatioCollection(hetAllelicCounts.getSampleName(), Collections.emptyList());     //create an empty collection with the appropriate name
        } else {
            readDenoisedCopyRatios();
            readAndFilterAllelicCounts();
            crafSegments = new MultidimensionalKernelSegmenter(denoisedCopyRatios, hetAllelicCounts)
                    .findSegmentation(maxNumSegmentsPerChromosome,
                            kernelVarianceCopyRatio, kernelVarianceAlleleFraction, 0.1, kernelApproximationDimension,
                            ImmutableSet.copyOf(windowSizes).asList(),
                            1., 1.);
        }
        writeSegments(crafSegments, CRAF_SEGMENTS_FILE_SUFFIX);

        logger.info("Modeling available denoised copy ratios and heterozygous allelic counts...");
        //initial MCMC model fitting performed by CRAFModeller constructor
        final AlleleFractionPrior alleleFractionPrior = new AlleleFractionPrior(minorAlleleFractionPriorAlpha);
        final CRAFModeller modeller = new CRAFModeller(
                crafSegments, denoisedCopyRatios, hetAllelicCounts, alleleFractionPrior,
                numSamplesCopyRatio, numBurnInCopyRatio,
                numSamplesAlleleFraction, numBurnInAlleleFraction);

        //write initial segments and parameters to file
        writeModeledSegmentsAndParameterFiles(modeller, BEGIN_FIT_FILE_TAG);

        //segmentation smoothing
        modeller.smoothSegments(
                maxNumSmoothingIterations, numSmoothingIterationsPerFit,
                smoothingCredibleIntervalThresholdCopyRatio, smoothingCredibleIntervalThresholdAlleleFraction);

        //write final segments and parameters to file
        writeModeledSegmentsAndParameterFiles(modeller, FINAL_FIT_FILE_TAG);

        logger.info("SUCCESS: ModelSegments run complete.");

        return "SUCCESS";
    }

    private void validateArguments() {
        Utils.nonNull(outputPrefix);
        Utils.validateArg(!(inputDenoisedCopyRatiosFile == null && inputAllelicCountsFile == null),
                "Must provide at least a denoised copy-ratio profile file or an allelic-counts file.");
        Utils.validateArg(!(inputAllelicCountsFile == null && inputNormalAllelicCountsFile != null),
                "Must provide an allelic-counts file for the case sample to run in matched-normal mode.");
        if (inputDenoisedCopyRatiosFile != null) {
            IOUtils.canReadFile(inputDenoisedCopyRatiosFile);
        }
        if (inputAllelicCountsFile != null) {
            IOUtils.canReadFile(inputAllelicCountsFile);
        }
        if (inputNormalAllelicCountsFile != null) {
            IOUtils.canReadFile(inputNormalAllelicCountsFile);
        }
    }

    private void readDenoisedCopyRatios() {
        logger.info(String.format("Reading denoised copy-ratio profile file (%s)...", inputDenoisedCopyRatiosFile));
        denoisedCopyRatios = new CopyRatioCollection(inputDenoisedCopyRatiosFile);
    }

    private CopyRatioSegmentCollection performCopyRatioSegmentation() {
        logger.info("Starting segmentation of denoised copy ratios...");
        final int maxNumChangepointsPerChromosome = maxNumSegmentsPerChromosome - 1;
        return new CopyRatioKernelSegmenter(denoisedCopyRatios)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVarianceCopyRatio, kernelApproximationDimension,
                        ImmutableSet.copyOf(windowSizes).asList(),
                        numChangepointsPenaltyFactorCopyRatio, numChangepointsPenaltyFactorCopyRatio);
    }

    private void readAndFilterAllelicCounts() {
        //read in case sample
        logger.info(String.format("Reading allelic-counts file (%s)...", inputAllelicCountsFile));
        final AllelicCountCollection unfilteredAllelicCounts = new AllelicCountCollection(inputAllelicCountsFile);
        final String sampleName = unfilteredAllelicCounts.getSampleName();

        //filter on total count in case sample
        logger.info(String.format("Filtering allelic counts with total count less than %d...", minTotalAlleleCount));
        AllelicCountCollection filteredAllelicCounts = new AllelicCountCollection(
                sampleName,
                unfilteredAllelicCounts.getRecords().stream()
                        .filter(ac -> ac.getTotalReadCount() >= minTotalAlleleCount)
                        .collect(Collectors.toList()));
        logger.info(String.format("Retained %d / %d sites after filtering on total count...",
                filteredAllelicCounts.getRecords().size(), unfilteredAllelicCounts.getRecords().size()));

        //filter on overlap with copy-ratio intervals, if available
        if (denoisedCopyRatios != null) {
            logger.info("Filtering allelic-count sites not overlapping with copy-ratio intervals...");
            final OverlapDetector<CopyRatio> copyRatioOverlapDetector = denoisedCopyRatios.getOverlapDetector();
            filteredAllelicCounts = new AllelicCountCollection(
                    sampleName,
                    filteredAllelicCounts.getRecords().stream()
                            .filter(copyRatioOverlapDetector::overlapsAny)
                            .collect(Collectors.toList()));
            logger.info(String.format("Retained %d / %d sites after filtering on overlap with copy-ratio intervals...",
                    filteredAllelicCounts.getRecords().size(), unfilteredAllelicCounts.getRecords().size()));
        }

        if (inputNormalAllelicCountsFile == null) {
            //filter on homozygosity in case sample
            logger.info("No matched normal was provided, not running in matched-normal mode...");
            logger.info("Performing binomial testing and filtering homozygous allelic counts...");
            hetAllelicCounts = new AllelicCountCollection(
                    unfilteredAllelicCounts.getSampleName(),
                    filteredAllelicCounts.getRecords().stream()
                            .filter(ac -> new BinomialTest().binomialTest(
                                    ac.getTotalReadCount(),
                                    Math.min(ac.getAltReadCount(), ac.getRefReadCount()),
                                    genotypingBaseErrorRate,
                                    AlternativeHypothesis.TWO_SIDED) < genotypingPValueThreshold + EPSILON)
                            .collect(Collectors.toList()));
            final File hetAllelicCountsFile = new File(outputDir, outputPrefix + HET_ALLELIC_COUNTS_FILE_SUFFIX);
            hetAllelicCounts.write(hetAllelicCountsFile);
            logger.info(String.format("Retained %d / %d sites after testing for heterozygosity...",
                    hetAllelicCounts.getRecords().size(), unfilteredAllelicCounts.getRecords().size()));
            logger.info(String.format("Heterozygous allelic counts written to %s.", hetAllelicCountsFile));
        } else {
            //read in matched normal
            logger.info("Matched normal was provided, running in matched-normal mode...");
            logger.info("Performing binomial testing and filtering homozygous allelic counts in matched normal...");
            final AllelicCountCollection unfilteredNormalAllelicCounts = new AllelicCountCollection(inputNormalAllelicCountsFile);
            if (!unfilteredNormalAllelicCounts.getIntervals().equals(unfilteredAllelicCounts.getIntervals())) {
                throw new UserException.BadInput("Allelic-count sites in case sample and matched normal do not match. " +
                        "Run CollectAllelicCounts using the same interval list of sites for both samples.");
            }
            final String normalSampleName = unfilteredNormalAllelicCounts.getSampleName();

            //filter on total count in matched normal
            logger.info(String.format("Filtering allelic counts in matched normal with total count less than %d...", minTotalAlleleCount));
            final AllelicCountCollection filteredNormalAllelicCounts = new AllelicCountCollection(
                    normalSampleName,
                    unfilteredNormalAllelicCounts.getRecords().stream()
                            .filter(ac -> ac.getTotalReadCount() >= minTotalAlleleCount)
                            .collect(Collectors.toList()));
            logger.info(String.format("Retained %d / %d sites in matched normal after filtering on total count...",
                    filteredNormalAllelicCounts.getRecords().size(), unfilteredNormalAllelicCounts.getRecords().size()));

            //filter on homozygosity in matched normal
            final AllelicCountCollection hetNormalAllelicCounts = new AllelicCountCollection(
                    normalSampleName,
                    filteredNormalAllelicCounts.getRecords().stream()
                            .filter(ac -> new BinomialTest().binomialTest(
                                    ac.getTotalReadCount(),
                                    Math.min(ac.getAltReadCount(), ac.getRefReadCount()),
                                    genotypingBaseErrorRate,
                                    AlternativeHypothesis.TWO_SIDED) < genotypingPValueThreshold + EPSILON)
                            .collect(Collectors.toList()));
            final File hetNormalAllelicCountsFile = new File(outputDir, outputPrefix + "." + normalSampleName + HET_ALLELIC_COUNTS_FILE_SUFFIX);
            hetNormalAllelicCounts.write(hetNormalAllelicCountsFile);
            logger.info(String.format("Retained %d / %d sites in matched normal after testing for heterozygosity...",
                    hetNormalAllelicCounts.getRecords().size(), unfilteredNormalAllelicCounts.getRecords().size()));
            logger.info(String.format("Heterozygous allelic counts for matched normal written to %s.", hetNormalAllelicCountsFile));

            //retrieve sites in case sample
            logger.info("Retrieving allelic counts at these sites in case sample...");
            final Set<SimpleInterval> hetNormalAllelicCountSites = new HashSet<>(hetNormalAllelicCounts.getIntervals());
            hetAllelicCounts = new AllelicCountCollection(
                    sampleName,
                    filteredAllelicCounts.getRecords().stream()
                            .filter(ac -> hetNormalAllelicCountSites.contains(ac.getInterval()))
                            .collect(Collectors.toList()));
            final File hetAllelicCountsFile = new File(outputDir, outputPrefix + HET_ALLELIC_COUNTS_FILE_SUFFIX);
            hetAllelicCounts.write(hetAllelicCountsFile);
            logger.info(String.format("Allelic counts for case sample at heterozygous sites in matched normal written to %s.", hetAllelicCountsFile));
        }
    }

    private AlleleFractionSegmentCollection performAlleleFractionSegmentation() {
        logger.info("Starting segmentation of heterozygous allelic counts...");
        final int maxNumChangepointsPerChromosome = maxNumSegmentsPerChromosome - 1;
        return new AlleleFractionKernelSegmenter(hetAllelicCounts)
                .findSegmentation(maxNumChangepointsPerChromosome, kernelVarianceAlleleFraction, kernelApproximationDimension,
                        ImmutableSet.copyOf(windowSizes).asList(),
                        numChangepointsPenaltyFactorAlleleFraction, numChangepointsPenaltyFactorAlleleFraction);
    }

    private void writeModeledSegmentsAndParameterFiles(final CRAFModeller modeller,
                                                       final String fileTag) {
        final ModeledSegmentCollection modeledSegments = modeller.getModeledSegments();
        writeSegments(modeledSegments, fileTag + SEGMENTS_FILE_SUFFIX);
        final File copyRatioParameterFile = new File(outputDir, outputPrefix + fileTag + CR_PARAMETER_FILE_SUFFIX);
        final File alleleFractionParameterFile = new File(outputDir, outputPrefix + fileTag + AF_PARAMETER_FILE_SUFFIX);
        modeller.writeModelParameterFiles(copyRatioParameterFile, alleleFractionParameterFile);
    }

    private void writeSegments(final LocatableCollection<?> segments,
                               final String fileSuffix) {
        final File segmentsFile = new File(outputDir, outputPrefix + fileSuffix);
        segments.write(segmentsFile);
        logger.info(String.format("Segments written to %s", segmentsFile));
    }
}