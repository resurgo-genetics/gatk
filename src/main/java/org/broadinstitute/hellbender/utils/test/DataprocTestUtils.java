package org.broadinstitute.hellbender.utils.test;

import org.broadinstitute.hellbender.utils.runtime.ProcessController;

import java.util.ArrayList;
import java.util.List;
import java.util.UUID;

/**
 * utilities to help run tests on a Dataproc cluster
 */
public final class DataprocTestUtils {

    private DataprocTestUtils(){}

    /**
     * create a minimal dataproc cluster which will shutdown automatically after 10 minutes of idling or 30 minutes
     * @return the name of the cluster
     */
    public static String createTestCluster() {
        final String clusterName = "gatk-test-" + UUID.randomUUID();
        final String[] command = new String[]{
                "gcloud", "beta", "dataproc", "clusters", "create",
                "--max-idle", "10m",
                "--max-age", "30m",
                "--num-workers", "2",
                "--master-machine-type", "n1-highmem-2",
                "--worker-machine-type", "n1-highmem-2",
                clusterName
        };
        BaseTest.runProcess(ProcessController.getThreadLocal(), command, "Couldn't create dataproc cluster");
        return clusterName;
    }

    /**
     * run a gatk spark tool on a dataproc cluster and throw if it fails
     * @param tool name of the tool
     * @param args arguments to the tool
     * @param clusterName name of the cluster
     */
    public static void launchGatkTool(final String tool, final List<String> args, final String clusterName){
            final List<String> command = new ArrayList<>();
            command.add("gatk-launch");
            command.add(tool);
            command.addAll(args);
            command.add("--");
            command.add("--sparkRunner"); command.add("GCS");
            command.add("--cluster"); command.add(clusterName);
        final String[] commandArray = command.toArray(new String[command.size()]);
        BaseTest.runProcess(ProcessController.getThreadLocal(), commandArray, "gatk run on dataproc failed");
    }
}
