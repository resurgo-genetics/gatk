package org.broadinstitute.hellbender.engine.spark;

import htsjdk.samtools.ValidationStringency;
import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.tools.spark.pipelines.PrintReadsSpark;
import org.broadinstitute.hellbender.utils.gcs.BucketUtils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.test.ArgumentsBuilder;
import org.broadinstitute.hellbender.utils.test.DataprocTestUtils;
import org.broadinstitute.hellbender.utils.test.IntegrationTestSpec;
import org.broadinstitute.hellbender.utils.test.SamAssertionUtils;
import org.testng.Assert;
import org.testng.annotations.BeforeClass;
import org.testng.annotations.DataProvider;
import org.testng.annotations.Test;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.Collections;

/**
 * tests that actuall run spark tools on dataproc
 */
public class DataprocIntegrationTest extends CommandLineProgramTest{
    private String clusterName;

    @BeforeClass(groups = "cloud")
    public void startCluster(){
        clusterName = DataprocTestUtils.createTestCluster();
    }

    @DataProvider
    public Object[][] getCloudPaths(){
        return new Object[][]{
                {"org/broadinstitute/hellbender/engine/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.10000000-10000020.with.unmapped.bam"},
                {"large/CEUTrio.HiSeq.WGS.b37.NA12878.20.21.bam"}
        };
    }

    @Test(dataProvider = "getCloudPaths", groups = "cloud")
    public void printReadSparkOnDataproc(String input) throws IOException {
        final String gcsInputPath = getGCPTestInputPath() + input;
        final String outputPath = BucketUtils.getTempFilePath(getGCPTestStaging(), ".bam");

        final ArgumentsBuilder argBuilder = new ArgumentsBuilder();
        argBuilder.addArgument("input", gcsInputPath)
                .addArgument("output", outputPath)
                .addArgument("bamPartitionSize", String.valueOf(10*1024*1024));
        DataprocTestUtils.launchGatkTool(PrintReadsSpark.class.getSimpleName(), argBuilder.getArgsList(), clusterName);

        final File expected = createTempFile("expected", ".bam");
        final File actual = createTempFile("actual",".bam");
        Assert.assertTrue(expected.delete());
        Assert.assertTrue(actual.delete());
        Files.copy(IOUtils.getPath(gcsInputPath), expected.toPath());
        Files.copy(IOUtils.getPath(outputPath), actual.toPath());
        IntegrationTestSpec.assertMatchingFiles(Collections.singletonList(actual), Collections.singletonList(expected.toString()), true, ValidationStringency.LENIENT);
    }
}
