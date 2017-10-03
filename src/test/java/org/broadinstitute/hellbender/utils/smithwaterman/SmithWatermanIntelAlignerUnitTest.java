package org.broadinstitute.hellbender.utils.smithwaterman;

public class SmithWatermanIntelAlignerUnitTest extends SmithWatermanAlignerAbstractUnitTest {

    /*
    *Test the Intel Aligner Native Interface
    */

    @Override
    protected SmithWatermanIntelAligner getAligner() {
        return new SmithWatermanIntelAligner();
    }


}
