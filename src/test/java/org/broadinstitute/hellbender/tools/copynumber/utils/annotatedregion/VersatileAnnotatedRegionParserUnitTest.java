package org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion;

import com.google.common.collect.Lists;
import org.broadinstitute.hellbender.utils.test.BaseTest;
import org.testng.Assert;
import org.testng.annotations.Test;
import org.testng.collections.Sets;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Set;

public class VersatileAnnotatedRegionParserUnitTest extends BaseTest {
    private final String TEST_FILE = publicTestDir + "org/broadinstitute/hellbender/tools/coveragemodel/learning_combined_copy_number.tsv";

    @Test
    public void basicTest() throws IOException {
        final VersatileAnnotatedRegionParser parser = new VersatileAnnotatedRegionParser();
        final Set<String> headersOfInterest = Sets.newHashSet(Lists.newArrayList("name", "learning_SAMPLE_0"));
        final List<SimpleAnnotatedGenomicRegion> simpleAnnotatedGenomicRegions =
                parser.readAnnotatedRegions(new File(TEST_FILE), headersOfInterest);

        Assert.assertEquals(simpleAnnotatedGenomicRegions.size(), 518);
        Assert.assertTrue(simpleAnnotatedGenomicRegions.stream()
                .mapToInt(s -> s.getAnnotations().entrySet().size())
                .allMatch(i -> i == headersOfInterest.size()));
        Assert.assertTrue(simpleAnnotatedGenomicRegions.stream().allMatch(s -> s.getAnnotations().keySet().containsAll(headersOfInterest)));

        // Grab a few at random and test values

    }
}
