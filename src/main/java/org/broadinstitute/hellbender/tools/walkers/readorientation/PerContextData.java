package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.broadinstitute.hellbender.engine.AlignmentContext;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.tools.walkers.readorientation.ContextDependentArtifactFilterEngine.State;

import java.util.ArrayList;
import java.util.List;

/**
 * Created by tsato on 8/7/17.
 */

/**
 * FIXME: perhaps a more natural approach is to organize data such that at index n we get all of the data for example n
 * Write now we have to get each field and then fastforward to nth example, which is probably wrong.
 */
public class PerContextData {

    private final String referenceContext;
    private final List<Nucleotide> altAlleles;
    private final List<String> positions; // for debugging

    // the number of data points (i.e. loci) with this 3-mer in the reference context
    private int numAltExamples = 0;
    private int numRefExamples = 0;

    private List<Integer> depths;
    private List<int[]> baseCounts;
    private List<int[]> f1r2Counts;

    final static int MAX_COVERAGE = 1000;
    final int[] refsiteCoverageHistogram = new int[MAX_COVERAGE]; // index 999 contains the count of sites with

    public PerContextData(final String referenceContext){
        this.referenceContext = referenceContext;

        // the list of most frequenctly observed alt allele at each site
        altAlleles = new ArrayList<>(ContextDependentArtifactFilter.DEFAULT_INITIAL_LIST_SIZE);
        positions = new ArrayList<>(ContextDependentArtifactFilter.DEFAULT_INITIAL_LIST_SIZE);

        depths = new ArrayList<>(ContextDependentArtifactFilter.DEFAULT_INITIAL_LIST_SIZE);
        baseCounts = new ArrayList<>(ContextDependentArtifactFilter.DEFAULT_INITIAL_LIST_SIZE);
        f1r2Counts = new ArrayList<>(ContextDependentArtifactFilter.DEFAULT_INITIAL_LIST_SIZE);
    }

    public void addAltExample(final int depth, final int[] altDepth, final int[] altF1R2Depth, final Nucleotide altAllele, final AlignmentContext alignmentContext){
        Utils.validateArg(MathUtils.sum(altDepth) >= MathUtils.sum(altF1R2Depth), "alt depths must be greater than or equal to alt F1R2 depth");
        depths.add(depth);
        baseCounts.add(altDepth);
        f1r2Counts.add(altF1R2Depth);
        altAlleles.add(altAllele); // TODO: do we need this?
        positions.add(String.format("%s:%s-%s", alignmentContext.getContig(), alignmentContext.getStart(), alignmentContext.getEnd()));

        numAltExamples++;
    }

    public void addRefExample(final int depth){
        // the last bin contains all sites that have coverage > MAX_COVERAGE
        final int index = depth < MAX_COVERAGE ? depth : MAX_COVERAGE - 1;
        refsiteCoverageHistogram[index] += 1;
        numRefExamples++;
    }

    // getters
    public String getReferenceContext(){ return referenceContext; }

    public List<Nucleotide> getAltAlleles(){ return altAlleles; }

    public int getNumExamples(){ return numRefExamples + numAltExamples; }

    public int getNumAltExamples(){ return numAltExamples; }

    public int getNumRefExamples(){ return numRefExamples; }

    public List<Integer> getDepths(){ return depths; }

    public List<int[]> getBaseCounts(){ return baseCounts; }

    public List<int[]> getF1r2Counts(){ return f1r2Counts; }

    public int[] getRefsiteCoverageHistogram() { return refsiteCoverageHistogram; }

    // TODO: write tests
    // return K-dimensional array of alt counts m_nk, returns a double[], not int[], because the downstream
    // API (MathArrays.ebeMultiply()) requires an array of doubles
    public double[] getAltCountsNicelyArranged(final int n){
        final Nucleotide altAllele = altAlleles.get(n);
        final double[] arrangedAltCounts = new double[State.values().length];
        final int[] baseCounts_n = baseCounts.get(n);
        final double altCountOfMostLikelyAllele = baseCounts_n[altAllele.ordinal()];
        final Nucleotide refAllele = Nucleotide.valueOf(referenceContext.substring(1, 2));
        final List<State> impossibleStates = State.getImpossibleStates(refAllele);

        for (State z : State.getBalancedStates()){
            arrangedAltCounts[z.ordinal()] = altCountOfMostLikelyAllele;
        }

        for (State z : State.artifactStates){
            if (impossibleStates.contains(z)){
                continue;
            }

            arrangedAltCounts[z.ordinal()] = baseCounts_n[State.getAltAlleleOfTransition(z).ordinal()];
        }
        return arrangedAltCounts;
    }

    // TODO: write tests
    // return K-dimensional array of alt f1r2 counts x_nk, similar to {@code getAltCountsNicelyArranged}
    public double[] getAltF1R2CountsNicelyArranged(final int n){
        final Nucleotide altAllele = altAlleles.get(n);
        final Nucleotide refAllele = Nucleotide.valueOf(referenceContext.substring(1, 2));
        final double[] arrangedAltF1R2Counts = new double[State.values().length];
        final int[] f1r2Counts_n = f1r2Counts.get(n);
        final double altF1R2CountOfMostLikelyAllele = f1r2Counts_n[altAllele.ordinal()];
        final List<State> impossibleStates = State.getImpossibleStates(refAllele);

        for (State z : State.getBalancedStates()){
            arrangedAltF1R2Counts[z.ordinal()] = altF1R2CountOfMostLikelyAllele;
        }

        for (State z : State.artifactStates){
            if (impossibleStates.contains(z)){
                continue;
            }

            arrangedAltF1R2Counts[z.ordinal()] = f1r2Counts_n[State.getAltAlleleOfTransition(z).ordinal()];
        }
        return arrangedAltF1R2Counts;
    }


}
