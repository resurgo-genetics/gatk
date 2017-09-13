package org.broadinstitute.hellbender.tools.walkers.readorientation;

import htsjdk.variant.variantcontext.Allele;
import org.apache.logging.log4j.LogManager;
import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.pileup.ReadPileup;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.broadinstitute.hellbender.tools.walkers.readorientation.ContextDependentArtifactFilterEngine.States;

import java.util.Arrays;

/**
 * Created by tsato on 8/14/17.
 */
public class ContextDependentArtifactFilterEngineUnitTest {
    private final byte A = "A".getBytes()[0];
    private final byte C = "C".getBytes()[0];
    private final byte G = "G".getBytes()[0];
    private final byte T = "T".getBytes()[0];

    private static final double EPSILON = 1e-4;
    /**
     * Create a test case: 100 sites for a single context AGT. G -> T transition observed in NUM_F1R2_EXAMPLES sites
     * All of the alt sites have 100% Alt F1R2 reads. The rest of the examples are hom ref.
     */
    @Test
    public void testSimpleCase() {
        final int NUM_EXAMPLES = 100;
        final int NUM_F1R2_EXAMPLES = 20;

        final int DEPTH = 100;
        final short STANDARD_ALT_DEPTH = 20;
        final short BALANCED_ALT_F1R2_DEPTH = STANDARD_ALT_DEPTH/2;
        final short ENTIRELY_F1R2 = STANDARD_ALT_DEPTH;

        final String refContext = "AGT";

        final PerContextData data = new PerContextData(refContext);
        for (int n = 0; n < NUM_EXAMPLES; n++){
            final SimpleInterval locus = new SimpleInterval("20", 1_000_000 + n, 1_000_000 + n);
            final AlignmentContext alignmentContext = new AlignmentContext(locus, new ReadPileup(locus));

            // assume 20 in 100 examples is alt, and 20% allele fraction, and alt reads are entirely F1R2
            if (n < NUM_F1R2_EXAMPLES){
                final short altDepth = STANDARD_ALT_DEPTH;
                final Nucleotide allele = Nucleotide.T;
                final short altF1R2Depth = ENTIRELY_F1R2;

                data.addAltExample(DEPTH, altDepth, altF1R2Depth, allele, alignmentContext);
            } else {
                data.addRefExample(DEPTH);
            }
        }

        ContextDependentArtifactFilterEngine engine = new ContextDependentArtifactFilterEngine(data);
        Hyperparameters hyperparameters = engine.runEMAlgorithm(LogManager.getLogger(this.getClass()));

        Assert.assertEquals(engine.effectiveCounts[States.F1R2.ordinal()], (double) NUM_F1R2_EXAMPLES, EPSILON);
        // START HERE, we get something like 78 - but we want 80. For such a simple example, we gotta get this right
        Assert.assertEquals(engine.effectiveCounts[States.BALANCED_HOM_REF.ordinal()], (double) NUM_EXAMPLES - NUM_F1R2_EXAMPLES, EPSILON);

        // Given that the allele is G, the ref allele, p(z = hom ref|a = G) must equal 1.0
        Assert.assertEquals(hyperparameters.getPi()[BaseUtils.simpleBaseToBaseIndex(G)][States.BALANCED_HOM_REF.ordinal()], 1.0, EPSILON);

        // Given that the allele is T, all of the examples were F1R2 artifacts, so we want p(z = F1R2|a = T) = 1.0 and
        // p(z = F2R1|a = T) = 0.0
        Assert.assertEquals(hyperparameters.getPi()[BaseUtils.simpleBaseToBaseIndex(T)][States.F1R2.ordinal()], 1.0);
        Assert.assertEquals(hyperparameters.getPi()[BaseUtils.simpleBaseToBaseIndex(T)][States.F2R1.ordinal()], 0.0);

        // We expect the model to learn the correct allele fraction given z = F1R2
        Assert.assertEquals(hyperparameters.getF()[States.F1R2.ordinal()], (double) ENTIRELY_F1R2/DEPTH);
    }


    /**
     * Now test the case where not all of the transitions have orientation bias. And on the transitions that do,
     * not all of the sites have orientation bias. Still assumes single context.
     */
    @Test
    public void testMoreComplicatedCase() {
        final String refContext = "AGT";

        final int NUM_EXAMPLES = 1000;
        final int NUM_F1R2_EXAMPLES_G_TO_T = 20;
        final int NUM_TOTAL_ALT_EXAMPLES = 100;

        final int DEPTH = 100;
        final short STANDARD_ALT_DEPTH = 20;
        final short BALANCED_ALT_F1R2_DEPTH = STANDARD_ALT_DEPTH/2;
        final short ENTIRELY_F1R2 = STANDARD_ALT_DEPTH;

        // first create the examples for the G -> T transitions, a fractin of which has read orientation bias
        final PerContextData data = new PerContextData(refContext);
        for (Nucleotide allele : Arrays.asList(Nucleotide.A, Nucleotide.C, Nucleotide.G, Nucleotide.T)) {
            for (int n = 0; n < NUM_EXAMPLES / 4; n++) {
                final SimpleInterval snpLocation = new SimpleInterval("20", 1_000_000 + n, 1_000_000 + n);
                final AlignmentContext alignmentContext = new AlignmentContext(snpLocation, new ReadPileup(snpLocation));

                if (n < NUM_TOTAL_ALT_EXAMPLES && allele != Nucleotide.G) {
                    final short altDepth = STANDARD_ALT_DEPTH;
                    final short altF1R2Depth = n < NUM_F1R2_EXAMPLES_G_TO_T && allele.equals(Nucleotide.T) ? ENTIRELY_F1R2 : BALANCED_ALT_F1R2_DEPTH;
                    data.addAltExample(DEPTH, altDepth, altF1R2Depth, allele, alignmentContext);
                } else {
                    // create hom ref examples
                    data.addRefExample(DEPTH);
                }




            }
        }

        // To consider: Can a site be both F1R2 and F2R1? What about the whole reverse complement thing?
        ContextDependentArtifactFilterEngine engine = new ContextDependentArtifactFilterEngine(data);
        Hyperparameters hyperparameters = engine.runEMAlgorithm(LogManager.getLogger(this.getClass()));

        Assert.assertEquals(engine.effectiveCounts[States.F1R2.ordinal()], (double) NUM_F1R2_EXAMPLES_G_TO_T, EPSILON);
        // temporariliy disable these two tests - first fix the error with pi
        Assert.assertEquals(engine.effectiveCounts[States.BALANCED_HOM_REF.ordinal()], (double) NUM_EXAMPLES - 3*NUM_TOTAL_ALT_EXAMPLES, EPSILON);
        Assert.assertEquals(engine.effectiveCounts[States.BALANCED_HET.ordinal()], (double) 3*NUM_TOTAL_ALT_EXAMPLES - NUM_F1R2_EXAMPLES_G_TO_T, EPSILON);

        // Given that the allele is G, or the ref allele, p(z = hom ref|a = G) must equal 1.0
        Assert.assertEquals(hyperparameters.getPi()[BaseUtils.simpleBaseToBaseIndex(G)][States.BALANCED_HOM_REF.ordinal()], 1.0, EPSILON);

        // Given that the allele is T, only fraction of the samples had the F1R2 artifacts
        Assert.assertEquals(hyperparameters.getPi()[BaseUtils.simpleBaseToBaseIndex(T)][States.F1R2.ordinal()], (double) NUM_F1R2_EXAMPLES_G_TO_T/NUM_TOTAL_ALT_EXAMPLES, EPSILON);
        Assert.assertEquals(hyperparameters.getPi()[BaseUtils.simpleBaseToBaseIndex(T)][States.BALANCED_HET.ordinal()], 1 - (double) NUM_F1R2_EXAMPLES_G_TO_T/NUM_TOTAL_ALT_EXAMPLES, EPSILON);

        // Given that the allele is A or C, all of the samples are balanced HET
        Assert.assertEquals(hyperparameters.getPi()[BaseUtils.simpleBaseToBaseIndex(A)][States.BALANCED_HET.ordinal()], 1.0, EPSILON);
        Assert.assertEquals(hyperparameters.getPi()[BaseUtils.simpleBaseToBaseIndex(C)][States.BALANCED_HET.ordinal()], 1.0, EPSILON);

        // We expect the model to learn the correct allele fractions
        Assert.assertEquals(hyperparameters.getF()[States.F1R2.ordinal()], (double) STANDARD_ALT_DEPTH/DEPTH, EPSILON);

        // TODO: if we learn theta, add a test here
    }

    /**
     * Test the case where the alt reads are heavily biased towards F1R2 or F2R1 but not entirely so
     */
    @Test
    public void testEvenMoreComplicated() {

    }

    /**
     * Test multiple contexts, not just the standard old
     */
    @Test
    public void testMultipleContexts() {

    }
}