package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.samtools.util.SequenceUtil;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.VariantProgramGroup;
import org.broadinstitute.hellbender.engine.*;
import org.broadinstitute.hellbender.engine.filters.ReadFilter;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.broadinstitute.hellbender.utils.read.ReadUtils;
import org.broadinstitute.hellbender.tools.walkers.readorientation.AltSiteRecord.AltSiteRecordTableWriter;
import org.broadinstitute.hellbender.tools.walkers.readorientation.RefSiteHistogram.RefSiteHistogramWriter;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Created by tsato on 7/26/17.
 */

@CommandLineProgramProperties(
        summary = "",
        oneLineSummary = "",
        programGroup = VariantProgramGroup.class
)

public class CollectDataForReadOrientationFilter extends LocusWalker {
    public static final String ALT_DATA_TABLE_SHORT_NAME = "alt_table";
    public static final String ALT_DATA_TABLE_LONG_NAME = "alt_data_table";

    public static final String REF_HISTOGRAM_TABLE_SHORT_NAME = "ref_table";
    public static final String REF_HISTOGRAM_TABLE_LONG_NAME = "ref_histogram_table";


    @Argument(fullName = "", shortName = "", doc = "exclude reads below this quality from pileup", optional = true)
    static int MINIMUM_MEDIAN_MQ = 20;

    @Argument(fullName = "", shortName = "", doc = "exclude bases below this quality from pileup", optional = true)
    static int MINIMUM_BASE_QUALITY = 10;

    @Argument(fullName = ALT_DATA_TABLE_LONG_NAME,
            shortName = ALT_DATA_TABLE_SHORT_NAME,
            doc = "a tab-separated table of data over alt sites")
    static File altDataTable = null;

    @Argument(fullName = REF_HISTOGRAM_TABLE_LONG_NAME,
            shortName = REF_HISTOGRAM_TABLE_SHORT_NAME,
            doc = "a tab-separated depth histogram over ref sites")
    static File refHistogramTable = null;

    private static RefSiteHistogram[] refSiteHistograms = new RefSiteHistogram[64];

    private static final List<String> ALL_3_MERS = SequenceUtil.generateAllKmers(3).stream().map(String::new).collect(Collectors.toList());

    public static final List<String> SOME_3_MERS = Arrays.asList("ACT", "GTT", "AAA", "CGT");

    AltSiteRecordTableWriter altTableWriter;

    RefSiteHistogramWriter refTableWriter;



    @Override
    public boolean requiresReference(){
        return true;
    }

    @Override
    public List<ReadFilter> getDefaultReadFilters() {
        return ReadUtils.makeStandardReadFilters();
    }

    @Override
    public void onTraversalStart() {
        for (final String refContext : ALL_3_MERS){
            refSiteHistograms[contextToIndex(refContext)] = new RefSiteHistogram(refContext);
        }

        // TODO: would it be faster to store the records in memory and write them all out in one go?
        // intentionally not use try-with-resources so that the writers stay open outside of the try block
        try {
            altTableWriter = new AltSiteRecordTableWriter(altDataTable);
        } catch (IOException e) {
            throw new UserException(String.format("Encountered an IO exception creating writers for %s or %s", altDataTable, refHistogramTable), e);
        }

    }

    @Override
    public void apply(final AlignmentContext alignmentContext, final ReferenceContext referenceContext, final FeatureContext featureContext){
        // referenceContext always comes with a window of a single base, so
        // manually expand the window and get the 3-mer for now.
        // TODO: implement getBasesInInterval() in referenceContext. Maybe simplify to getKmer(int k)?
        // TODO: this is still relevant (10/2). I shouldn't mess with the internal state of the ref context object
        referenceContext.setWindow(1, 1);
        final String refContext = new String(referenceContext.getBases());
        assert refContext.length() == 3 : "kmer must have length 3";
        if (refContext.contains("N")) {
            return;
        }

        if (refContext == null){
            logger.info(String.format("null reference found at interval %s, k-mer = %s",
                    referenceContext.getInterval().toString(), refContext));
            return;
        }

        final ReadPileup pileup = alignmentContext.filter(pe -> pe.getQual() > MINIMUM_BASE_QUALITY).getBasePileup();

        // FIXME; this is not ideal. AlignmentContext should come filtered and not reach here if it's empty
        if (pileup.size() == 0){
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

        final Nucleotide refBase = Nucleotide.valueOf(refContext.getBytes()[1]);

        // make a copy fo base counts, update the counts of ref to -infty. Now the argmax of the array gives us
        // the alt base.
        final int[] baseCountsCopy = Arrays.copyOf(baseCounts, baseCounts.length);
        baseCountsCopy[refBase.ordinal()] = Integer.MIN_VALUE;
        final int altBaseIndex = MathUtils.argmax(baseCountsCopy);
        final boolean referenceSite = baseCounts[altBaseIndex] == 0;

        /*** End heuristics ***/

        // if the site is ref, we simply update the coverage histogram
        if (referenceSite){
            refSiteHistograms[contextToIndex(refContext)].increment(depth);
            return;
        }

        // we have an alt site
        final Nucleotide altBase = Nucleotide.valueOf(BaseUtils.baseIndexToSimpleBase(altBaseIndex));

        // TODO: this may be refactored
        final int[] altF1R2Counts = new int[4];
        for (Nucleotide nucleotide : Arrays.asList(Nucleotide.A, Nucleotide.C, Nucleotide.G, Nucleotide.G)){
            altF1R2Counts[nucleotide.ordinal()] = pileup.getNumberOfElements(pe -> Nucleotide.valueOf(pe.getBase()) == nucleotide && ! ReadUtils.isF2R1(pe.getRead()));
        }

        try {
            altTableWriter.writeRecord(new AltSiteRecord(refContext, baseCounts, altF1R2Counts, depth, altBase));
        } catch (IOException e) {
            throw new UserException("Encountered an IO Exception writing to the alt data table", e);
        }

        return;
    }

    @Override
    public Object onTraversalSuccess() {
        RefSiteHistogram.writeRefSiteHistograms(Arrays.asList(refSiteHistograms), refHistogramTable);
        return "SUCCESS";
    }

    @Override
    public void closeTool() {
        if (altTableWriter != null) {
            try {
                altTableWriter.close();
            } catch (IOException e) {
                throw new UserException("Encountered an IO exception while closing the alt table writer", e);
            }
        }

        if (refTableWriter != null) {
            try {
                refTableWriter.close();
            } catch (IOException e) {
                throw new UserException("Encountered an IO exception while closing the ref table writer", e);
            }
        }
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
