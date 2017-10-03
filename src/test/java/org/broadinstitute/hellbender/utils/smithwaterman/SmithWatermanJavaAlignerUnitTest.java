package org.broadinstitute.hellbender.utils.smithwaterman;

import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.testng.Assert;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

public final class SmithWatermanJavaAlignerUnitTest extends SmithWatermanAlignerAbstractUnitTest {

    @Override
    protected SmithWatermanJavaAligner getAligner() {
        return SmithWatermanJavaAligner.getInstance();
    }




}
