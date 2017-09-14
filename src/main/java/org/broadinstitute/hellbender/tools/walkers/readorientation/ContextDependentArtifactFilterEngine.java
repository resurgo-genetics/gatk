package org.broadinstitute.hellbender.tools.walkers.readorientation;

import org.apache.commons.math3.distribution.BinomialDistribution;
import org.apache.commons.math3.util.MathArrays;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.utils.BaseUtils;
import org.broadinstitute.hellbender.utils.GATKProtectedMathUtils;
import org.broadinstitute.hellbender.utils.MathUtils;
import org.broadinstitute.hellbender.utils.Nucleotide;
import org.spark_project.guava.annotations.VisibleForTesting;

import java.util.Arrays;
import java.util.stream.IntStream;

/**
 * Created by tsato on 7/26/17.
 */
public class ContextDependentArtifactFilterEngine {
    // z \in { F1R2, F2R1, Balanced Hom Ref, Balanced Het, Balanced Hom Var }. Thus |z| = K = 5.
    static final int NUM_STATUSES = 5;

    static final String[] ALL_ALLELES = new String[] { "A", "C", "G", "T" };


    // A, C, G, or T, since we only look at SNP sites
    static final int NUM_ALLELES = ALL_ALLELES.length; // aka 4

    // When the increase in likelihood falls under this value, we call the algorithm converged
    static final double CONVERGENCE_THRESHOLD = 1e-3;

    // Regularizer (?) TODO: think this through
    static final double EPSILON = 1e-4;


    final String referenceContext;

    final Nucleotide refAllele;

    // Observed data
    final PerContextData data;

    // N by K matrix of (natural) log posterior probabilities of latent variable z
    // with the current estimates of hyperparameters pi, f, and theta, where N is the number of alt sites
    // under the context
    final double[][] log10AltResponsibilities;

    // experimental measure - eventually pick a log or linear version
    final double[][] altResponsibilities;

    // MAX_VALUE by K matrix of a cache of responsibilities of a ref site (i.e. m = 0, x = 0)
    // for ref sites with coverage 0, 1, ..., MAX_VALUE - 1.
    final double[][] refResponsibilities = new double[PerContextData.MAX_COVERAGE][NUM_STATUSES];

    final int numExamples;

    // A by K matrix of effective counts. Each column must add up to N_k
    /**
     *              State z
     *             _ _ _ _ _
     *         A |
     * alleles C |
     *         G |
     *         T |
     **/
    final double[][] effectiveCountsGivenAllele = new double[NUM_ALLELES][NUM_STATUSES];

    // K-dimensional vector of effective sample counts for each class of z, weighted by the the altResponsibilities. For a fixed k,
    // we sum up the counts over all alleles. N_k in the docs.
    // TODO: should be final
    @VisibleForTesting
    double[] effectiveCounts = new double[NUM_STATUSES];

    // A-dimensional vector of effective count of samples for each allele, where A = 4 = |{A, C, G, T}|, N_a in docs
    @VisibleForTesting
    double[] effectiveCountsOfAlleles = new double[NUM_ALLELES];

    /*** Hyperparameters of the model ***/

    // pi is the A by K matrix of weights (prior?) for the categorical variable z
    // For each allele a, the K-dimensional vector pi_a must be a valid probability distribution over z
    // In other words, each row must add up to 1
    final double[][] pi = new double[NUM_ALLELES][NUM_STATUSES];

    // K-dimensional vector of probabilities for the binomial m given z, which represents the number of alt reads at site n
    final double[] f = new double[NUM_STATUSES];

    // In case we observe no data - assume that the allele fraction given artifact is this value
    static final double DEFAULT_ARTIFACT_ALLELE_FRACTION = 0.3;

    // K-dimensional vector of probabilities for the binomial x given z, which represents the number of F1R2 alt reads at site n
    final double[] theta = new double[NUM_STATUSES];

    private int numIterations = 0;

    private static int MAX_ITERATIONS = 200;

    // one may plot the changes in L2 distance of parameters to make sure that EM is steadily moving towards the (local? global?) maximum
    // TODO: this being public is questionable
    public double[] l2distancesOfParameters = new double[MAX_ITERATIONS];


    public ContextDependentArtifactFilterEngine(final PerContextData data){
        this.data = data;
        numExamples = data.getNumExamples();
        log10AltResponsibilities = new double[numExamples][NUM_STATUSES];
        altResponsibilities = new double[numExamples][NUM_STATUSES];

        // initialize altResponsibilities in log 10 space
        final double initialLog10Probability = - Math.log10(NUM_STATUSES);
        for (int n = 0; n < data.getNumExamples(); n++ ) {
            Arrays.fill(log10AltResponsibilities[n], initialLog10Probability);
        }

        final double initialProbability = 1.0/NUM_STATUSES;
        for (int n = 0; n < data.getNumExamples(); n++ ) {
            Arrays.fill(altResponsibilities[n], initialProbability);
        }

        for (int n = 0; n < PerContextData.MAX_COVERAGE; n++ ) {
            Arrays.fill(refResponsibilities[n], initialProbability);
        }

        referenceContext = data.getReferenceContext();
        refAllele = Nucleotide.valueOf(referenceContext.substring(1,2));

        // populate alleles
        final String referenceBase = referenceContext.substring(1, 2);

        // we fix some of the parameters to entice the model to assign particular states to the indices into z
        // for instance, we fix the allele fraction parameter f for z = Balanced Het to be 0.5.
        f[States.BALANCED_HET.ordinal()] = 0.5;
        f[States.BALANCED_HOM_REF.ordinal()] = EPSILON;
        f[States.BALANCED_HOM_VAR.ordinal()] = 1 - EPSILON;

        // similarly, we may fix some of theta_z
        // TODO: or should we learn them?
        theta[States.F1R2.ordinal()] = 1 - EPSILON;
        theta[States.F2R1.ordinal()] = EPSILON;
        theta[States.BALANCED_HOM_REF.ordinal()] = 0.5;
        theta[States.BALANCED_HET.ordinal()] = 0.5;
        theta[States.BALANCED_HOM_VAR.ordinal()] = 0.5;
    }

    public Hyperparameters runEMAlgorithm(final Logger logger){
        boolean converged = false;
        double[][] oldPi = new double[NUM_ALLELES][NUM_STATUSES];

        for (int a = 0; a < NUM_ALLELES; a++){
             oldPi[a] = Arrays.copyOf(pi[a], NUM_STATUSES);
        }


        while (!converged && numIterations < MAX_ITERATIONS){
            // TODO: stylistic problems here, there's too much side-effect stuff
            takeMstep();

            // assert newLikelihood >= oldLikelihood : "M step must increase the likelihood";
            final double l2Distance = IntStream.range(0, NUM_ALLELES).mapToDouble(a -> MathArrays.distance(oldPi[a], pi[a])).sum();
            converged = l2Distance < CONVERGENCE_THRESHOLD;

            l2distancesOfParameters[numIterations] = l2Distance;

            for (int a = 0; a < NUM_ALLELES; a++){
                oldPi[a] = Arrays.copyOf(pi[a], NUM_STATUSES);
            }

            takeEstep();

            numIterations++;
        }

        logger.info(String.format("Context %s, EM converged in %d steps", referenceContext, numIterations));
        return new Hyperparameters(referenceContext, pi, f, theta);
    }

    // Given the current estimates of the parameters pi, f, and theta, compute the log10AltResponsibilities
    // gamma_nk = p(z_nk|a)
    private void takeEstep(){
        // all ref sites look the same: altDepth = 0, altF1R2Depth = 0
        // When we have no alt reads (i.e. m = 0), the binomial over the number of F1R2 is a deterministic;
        // namely, it's always 1, and log(1) = 0
        for (int depth = 0; depth < PerContextData.MAX_COVERAGE; depth++){
            final int m = depth; // another hack to use depth in a stream
            final double[] log10UnnormalizedResponsibilities = IntStream.range(0, NUM_STATUSES)
                    .mapToDouble(k -> Math.log(pi[refAllele.ordinal()][k]) * MathUtils.LOG10_OF_E + new BinomialDistribution(m, f[k]).logProbability(0) * MathUtils.LOG10_OF_E)
                    .toArray();
            refResponsibilities[depth] = MathUtils.normalizeFromLog10ToLinearSpace(log10UnnormalizedResponsibilities);
        }

        // we must compute the altResponsibilities of each of n alt sites \gamma_{nk}
        for (int example = 0; example < data.getNumAltExamples(); example++){
            final int n = example; // hack to work around the fact that java stream doesn't let you use a non-final variable

            final Nucleotide allele = data.getAlleles().get(n);
            final int depth = data.getDepths().get(n);
            final short altDepth = data.getAltDepths().get(n);
            final short altF1R2Depth = data.getAltF1R2Depths().get(n);
            assert altDepth > 0 : "somehow a ref site slipped in as alt";

            final int a = BaseUtils.simpleBaseToBaseIndex(allele.toBase());

            final double[] log10AlleleFractionTerms = IntStream.range(0, NUM_STATUSES)
                    .mapToDouble(k -> new BinomialDistribution(depth, f[k]).logProbability(altDepth) * MathUtils.LOG10_OF_E)
                    .toArray();

            final double[] logAltF1R2FractionTerms = IntStream.range(0, NUM_STATUSES)
                    .mapToDouble(k -> new BinomialDistribution(altDepth, theta[k]).logProbability(altF1R2Depth) * MathUtils.LOG10_OF_E)
                    .toArray();

            assert log10AlleleFractionTerms.length == NUM_STATUSES : "alleleFractionFactors must have length K";
            assert logAltF1R2FractionTerms.length == NUM_STATUSES : "altF1R2FractionFactors must have length K";

            // TODO: might want to do this in log space, watch out for underflow.
            double[] log10UnnormalizedResponsibilities = IntStream.range(0, NUM_STATUSES)
                    .mapToDouble(k -> Math.log10(pi[a][k]) + log10AlleleFractionTerms[k] + logAltF1R2FractionTerms[k])
                    .toArray();

            // do we need the log 10 responsibilities?
            log10AltResponsibilities[n] = MathUtils.normalizeLog10(log10UnnormalizedResponsibilities, true, false);

            // do we really need to normalize here?
            altResponsibilities[n] = MathUtils.normalizeFromLog10ToLinearSpace(log10UnnormalizedResponsibilities);

            assert Math.abs(MathUtils.sumLog10(log10AltResponsibilities[n]) - 1.0) < EPSILON :
                    String.format("log responsibility for %dth example added up to %f", n,  MathUtils.sumLog10(log10AltResponsibilities[n]));
            assert Math.abs(MathUtils.sum(altResponsibilities[n]) - 1.0) < EPSILON :
                    String.format("responsibility for %dth example added up to %f", n,  MathUtils.sumLog10(altResponsibilities[n]));
        }
    }

    // given the current posterior distributions (i.e. altResponsibilities) over z, compute the estimate for
    // the categorical weights (pi), allele fractions (f), and alt F1R2 fraction (theta) that maximizes the lower bound
    // for the marginal likelihood, which is equivalent to maximizing the expectation of the complete data likelihood
    // with respect to the posterior over z
    private void takeMstep(){
        /*** compute responsibility-based statistics based on the current log10AltResponsibilities ***/

        // reset the effectiveCountsGivenAllele array
        for (int a = 0; a < NUM_ALLELES; a++){
            Arrays.fill(effectiveCountsGivenAllele[a], 0.0);
        }

        // TODO: optimize
        for (int n = 0; n < data.getNumAltExamples(); n++) {
            final Nucleotide allele = data.getAlleles().get(n);
            final int a = BaseUtils.simpleBaseToBaseIndex(allele.toBase()); // allele.ordinal()?

            // Warning; raising the log altResponsibilities by 10 may be naive REWORD
            effectiveCountsGivenAllele[a] = MathArrays.ebeAdd(effectiveCountsGivenAllele[a],
                    Arrays.stream(log10AltResponsibilities[n]).map(logp -> Math.pow(10, logp)).toArray());
        }

        // we optimize the sum of responsibilities over ref allele by taking a dot product of the cached responsibilities
        // at depths [0, 1, 2, ...., MAX_COVERAGE] and the number of sites at those depths
        effectiveCountsGivenAllele[refAllele.ordinal()] = GATKProtectedMathUtils.sumArrayFunction(0, PerContextData.MAX_COVERAGE,
                d -> MathArrays.scale(data.getRefsiteCoverageHistogram()[d], refResponsibilities[d]));
        assert Math.abs(MathUtils.sum(effectiveCountsGivenAllele[refAllele.ordinal()]) - data.getNumRefExamples()) < EPSILON :
                String.format("Effective count given ref allele should equal the number of ref sites (%d) but got %.3f",
                        data.getNumRefExamples(),
                        MathUtils.sum(effectiveCountsGivenAllele[refAllele.ordinal()]));

        // For a fixed k, sum up the effective counts over all alleles. N_k in the docs.
        effectiveCounts = GATKProtectedMathUtils.sumArrayFunction(0, NUM_ALLELES, a -> effectiveCountsGivenAllele[a]);

        assert effectiveCounts.length == NUM_STATUSES : "effectiveCount must be a k-dimensional vector";
        assert Math.abs(MathUtils.sum(effectiveCounts) - numExamples) < EPSILON :
                String.format("effective counts must add up to number of examples %d but got %f", numExamples, MathUtils.sum(effectiveCounts));

        // TODO: we don't have a good way of adding up columns of a 2-dimensional array
        for (int a = 0; a < NUM_ALLELES; a++){
            effectiveCountsOfAlleles[a] = MathUtils.sum(effectiveCountsGivenAllele[a]);
        }

        assert Math.abs(MathUtils.sum(effectiveCountsOfAlleles) - numExamples) < EPSILON : "effectiveCountOfAlleles should add up to numExamples";

        // K-dimensional vector of sample means weighted by the log10AltResponsibilities, N_k \bar{m} in the docs
        // ref sites do not contribute to the sum as m_n = 0 for ref sites
        // TOOD: I should probabily not enter log space - where do we need it? Where do we multiply a lot of small numbers?
        final double[] weightedSampleMeanM = GATKProtectedMathUtils.sumArrayFunction(0, data.getNumAltExamples(),
                n -> MathArrays.scale((double) data.getAltDepths().get(n), altResponsibilities[n]));
        assert weightedSampleMeanM.length == NUM_STATUSES : "weightedSampleMeanM should have length K";

        // K-dimensional vector of sample means weighted by the log10AltResponsibilities, N_k \bar{x} in the docs
        // ref sites do not contribute to the sum as x_n = 0 for ref sites
        final double[] weightedSampleMeanX = GATKProtectedMathUtils.sumArrayFunction(0, data.getNumAltExamples(),
                n -> MathArrays.scale((double) data.getAltF1R2Depths().get(n), altResponsibilities[n]));

        // K-dimensional vector of mean read depths weighted by the log10AltResponsibilities, N_k \bar{R} in the docs
        final double[] weightedDepthRAlt = GATKProtectedMathUtils.sumArrayFunction(0, data.getNumAltExamples(),
                n -> MathArrays.scale((double) data.getDepths().get(n), altResponsibilities[n]));
        final double[] weightedDepthRRef = GATKProtectedMathUtils.sumArrayFunction(0, PerContextData.MAX_COVERAGE,
                d -> MathArrays.scale(data.refsiteCoverageHistogram[d] * d, refResponsibilities[d]));
        final double[] weightedDepthR = MathArrays.ebeAdd(weightedDepthRAlt, weightedDepthRRef);

        // TODO: should we learn theta?
        // theta = MathArrays.ebeDivide(weightedSampleMeanX, weightedSampleMeanM);

        // We update some allele fractions according to data and log10AltResponsibilities and keep others fixed
        final double[] updatedAlleleFractions = MathArrays.ebeDivide(weightedSampleMeanM, weightedDepthR);
        f[States.F1R2.ordinal()] = updatedAlleleFractions[States.F1R2.ordinal()] == 0 ? DEFAULT_ARTIFACT_ALLELE_FRACTION :
                updatedAlleleFractions[States.F1R2.ordinal()];
        f[States.F2R1.ordinal()] = updatedAlleleFractions[States.F2R1.ordinal()]  == 0 ? DEFAULT_ARTIFACT_ALLELE_FRACTION :
                updatedAlleleFractions[States.F2R1.ordinal()];

        // update pi
        for (int a = 0; a < NUM_ALLELES; a++){
            // N_a in the docs
            final double numExamplesWithThisAllele = effectiveCountsOfAlleles[a];
            // when there are no examples of this particular (context, allele) pair, we give pi a uniform distribution
            if (numExamplesWithThisAllele < EPSILON){
                Arrays.fill(pi[a], 1/NUM_STATUSES);
            } else {
                pi[a] = MathArrays.scale(1 / numExamplesWithThisAllele, effectiveCountsGivenAllele[a]);
            }

            // ensure that each row of the pi matrix is normalized
            final double sumProbabilities = Math.abs(MathUtils.sum(pi[a]));
            assert sumProbabilities - 1.0 < EPSILON : "pi[a] must add up to 1";
        }

        return;
    }

    private boolean checkLikelihoodHasConverged(final double oldLikelihood, final double newLikelihood){
        return Math.abs(newLikelihood - oldLikelihood) < CONVERGENCE_THRESHOLD;
    }

    private static boolean checkLikelihoodHasConverged(final int numIterations){
        return numIterations > 10;
    }

    public enum States {
        F1R2, // Orientation bias at the site, evidence for alt is predominantly F1R2 reads
        F2R1, // Orientation Bias at the site, evidence for alt is predominantly F1R1 reads
        BALANCED_HOM_REF, // No orientation bias, and the site is hom ref
        BALANCED_HET, // No orientation bias, and the site is het
        BALANCED_HOM_VAR // No orientation bias, and the site is hom var
    }

}
