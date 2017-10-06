package org.broadinstitute.hellbender.utils.smithwaterman;


import com.intel.gkl.smithwaterman.IntelSmithWaterman;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWAlignerNativeBinding;
import org.broadinstitute.hellbender.exceptions.UserException;


public final class SmithWatermanIntelAligner implements SmithWatermanAligner {

    private static final SmithWatermanIntelAligner ALIGNER = new SmithWatermanIntelAligner();

    /**
     * return the stateless singleton instance of SmithWatermanIntelAligner
     */
    public static SmithWatermanIntelAligner getInstance() {
        return ALIGNER;
    }

    private final SWAlignerNativeBinding SmithWaterman = new IntelSmithWaterman();

    /*
    * Generate SWAlignerWrapper instance
    */
    private final SWNativeAlignerWrapper alignerWrapper = new SWNativeAlignerWrapper(SmithWaterman);


    /**
     * Create a new SW native pairwise aligner
     */
    public SmithWatermanIntelAligner() {
        final boolean isSupported = SmithWaterman.load(null);
        if (!isSupported) {
            throw new UserException.HardwareFeatureException("Machine does not support AVX SmithWaterman.");
        }
    }

    /**
     * Aligns the alternate sequence to the reference sequence
     *
     * @param reference  ref sequence
     * @param alternate  alt sequence
     */
    @Override
    public SmithWatermanAlignment align(final byte[] reference, final byte[] alternate, final SWParameters parameters, final SWOverhangStrategy overhangStrategy) {
        return alignerWrapper.align(reference, alternate,parameters,overhangStrategy);
    }

    /**
     * Close the aligner
     */
    @Override
    public void close() {
        alignerWrapper.close();
    }
}
