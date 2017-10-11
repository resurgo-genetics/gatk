package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ReadUtils;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by tsato on 7/26/17.
 */

/***
 * This tools is the learning phase of the orientation filter.
 * Inference phase will likely feature variant context and what not.
 */
@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = VariantProgramGroup.class
)

public class ContextDependentArtifactFilter extends LocusWalker {
    @Argument(fullName = "", shortName = "", doc = "exclude reads below this quality from pileup", optional = true)
    static int MINIMUM_MEDIAN_MQ = 20;

    @Argument(fullName = "", shortName = "", doc = "exclude bases below this quality from pileup", optional = true)
    static int MINIMUM_BASE_QUALITY = 10;

    @Argument(fullName = "test", shortName = "test", doc = "", optional = true)
    static boolean test = false;

    @Argument(fullName = StandardArgumentDefinitions.OUTPUT_LONG_NAME, shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME, doc = "a tab-seprated table of hyperparameters")
    static File output = null;

    public static PerContextData[] contextDependentDataMap;

    public static final List<String> ALL_3_MERS = SequenceUtil.generateAllKmers(3).stream().map(String::new).collect(Collectors.toList());

    public static final List<String> SOME_3_MERS = Arrays.asList("ACT", "GTT", "AAA", "CGT");

    // 37M bases in the exome, and there's about 1 SNP every 1000 bases. If we assume all 64 contexts are equally likely
    // we get 37_000/64 = 57. Use 64 for the closest power of 2
    static final int DEFAULT_INITIAL_LIST_SIZE = 64;


    @Override
    public boolean requiresReference(){
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return ReadUtils.makeStandardReadFilters();
    }

    @Override
    public void onTraversalStart(){
        contextDependentDataMap = new PerContextData[64]; // 4^3 = 64

        for (final String refContext : ALL_3_MERS){
            contextDependentDataMap[contextToIndex(refContext)] = new PerContextData(refContext);
        }
    }

    @Override
    public void apply(final AlignmentContext alignmentContext, final ReferenceContext referenceContext, final FeatureContext featureContext){
        // referenceContext always comes withe window of single base, so
        // manually expand the window and get the 3-mer for now.
        // TODO: implement getBasesInInterval() in referenceContext. Maybe simplify to getKmer(int k)?
        // TODO: this is still relevant (10/2). I shouldn't mess with the internal state of the ref context object
        referenceContext.setWindow(1, 1);
        final String reference3mer = new String(referenceContext.getBases());
        assert reference3mer.length() == 3 : "kmer must have length 3";
        if (reference3mer.contains("N")) {
            return;
        }

        if (reference3mer == null){
            logger.info(String.format("null reference found at interval %s, k-mer = %s",
                    referenceContext.getInterval().toString(), reference3mer));
            return;
        }

        final ReadPileup pileup = alignmentContext.filter(pe -> pe.getQual() > MINIMUM_BASE_QUALITY).getBasePileup();

        // FIXME; this is not ideal. AlignmentContext should come filtered and not reach here if it's empty
        if (pileup.size() == 0){
            logger.info(String.format("Empty pileup at position %s:%d-%d",
                    alignmentContext.getContig(), alignmentContext.getStart(), alignmentContext.getEnd()));
            return;
        }

        /*** Start heuristics ***/

        // skip INDELs

        // skip MQ=0 loci
        List<Integer> mappingQualities = new ArrayList<>(pileup.size());

        // there is no shortcut or a standard API for converting an int[] to List<Integer> (we don't want List<int[]>)
        // so we must convert int[] to a List<Integer> with a for loop
        // Median in Apache commons takes in a double[], not int[], so it's not much of an improvement
        for (final int mq : pileup.getMappingQuals()) {
            mappingQualities.add(mq);
        }

        final double medianMQ = MathUtils.median(mappingQualities);

        if (medianMQ < MINIMUM_MEDIAN_MQ) {
            return;
        }

        final int[] baseCounts = pileup.getBaseCounts();

        // R in the docs
        final int depth = (int) MathUtils.sum(baseCounts);

        final Nucleotide refBase = Nucleotide.valueOf(reference3mer.getBytes()[1]);

        // make a copy fo base counts, update the counts of ref to -infty. Now the argmax of the array gives us
        // the alt base.
        final int[] baseCountsCopy = Arrays.copyOf(baseCounts, baseCounts.length);
        baseCountsCopy[refBase.ordinal()] = Integer.MIN_VALUE;
        final int altBaseIndex = MathUtils.argmax(baseCountsCopy);
        final boolean referenceSite = baseCounts[altBaseIndex] == 0;

        /*** End heuristics ***/

        // if the site is ref, we simply update the coverage histogram
        if (referenceSite){
            // TODO: accessing the map excessively. Use array, it'd be faster
            contextDependentDataMap[contextToIndex(reference3mer)].addRefExample(depth);
            return;
        }

        // we have an alt site
        final Nucleotide altBase = Nucleotide.valueOf(BaseUtils.baseIndexToSimpleBase(altBaseIndex));

        // m in the docs
        final short altDepth = (short) baseCounts[altBase.ordinal()];

        // x in the docs
        final short altF1R2Depth = (short) pileup.getNumberOfElements(pe -> Nucleotide.valueOf(pe.getBase()) == altBase && ! ReadUtils.isF2R1(pe.getRead()));

        // TODO: this may be refactored
        final int[] altF1R2Counts = new int[4];
        for (Nucleotide nucleotide : Arrays.asList(Nucleotide.A, Nucleotide.C, Nucleotide.G, Nucleotide.G)){
            altF1R2Counts[nucleotide.ordinal()] = pileup.getNumberOfElements(pe -> Nucleotide.valueOf(pe.getBase()) == nucleotide && ! ReadUtils.isF2R1(pe.getRead()));
        }

        contextDependentDataMap[contextToIndex(reference3mer)].addAltExample(depth, baseCounts, altF1R2Counts, altBase, alignmentContext);
        return;
    }

    @Override
    public Object onTraversalSuccess() {
        List<Hyperparameters> hyperparameterEstimates = new ArrayList<>();

        // remember we run EM separately for each of 4^3 = 64 ref contexts
        for (final String refContext : test ? SOME_3_MERS : ALL_3_MERS){
            final PerContextData data = contextDependentDataMap[contextToIndex(refContext)];
            if (data.getNumAltExamples() == 0){
                // without alt examples there's no point in reporting any probabilities
                continue;
            }

            final ContextDependentArtifactFilterEngine engine = new ContextDependentArtifactFilterEngine(data);
            final Hyperparameters hyperparameters = engine.runEMAlgorithm(logger);
            hyperparameterEstimates.add(hyperparameters);
        }

        Hyperparameters.writeHyperparameters(hyperparameterEstimates, output);

        return "SUCCESS";
    }

    /***
     * Maps a reference 3-mer to an array index using the quaternary (base-4) numeral system
     * Example: AGT is represented as 023 in quaternary, which in decimal is  2*4^1 + 3*4^0 = 8+3 = 11
     *
     * @param reference3mer
     * @return
     */
    protected int contextToIndex(String reference3mer){
        Utils.validateArg(reference3mer.matches("[ACGT]{3}"),
                "input must be a string of length 3 and comprise of A, C, G,and T");
        final int digit2 = Nucleotide.valueOf(reference3mer.substring(0, 1)).ordinal();
        final int digit1 = Nucleotide.valueOf(reference3mer.substring(1, 2)).ordinal();
        final int digit0 = Nucleotide.valueOf(reference3mer.substring(2, 3)).ordinal();

        return digit2 * 16 + digit1 * 4 + digit0;
    }
}
