package org.broadinstitute.hellbender.tools.copynumber;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.formats.CopyNumberStandardArgument;
import org.junit.Test;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public class ModelSegmentsIntegrationTest extends CommandLineProgramTest {
    @Test
    public void testWGSChr20Chr21() {
        final String[] arguments = {
                "-" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-01A-01D-1931-08-chr20-chr21.denoisedCR.tsv",
                "-" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-01A-01D-1931-08-chr20-chr21.allelicCounts.tsv",
                "-" + CopyNumberStandardArgument.OUTPUT_PREFIX_SHORT_NAME, "TCGA-05-4389-01A-01D-1931-08-chr20-chr21",
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "/home/slee/working/gatk",
                "-" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testWGSChr20Chr21Normal() {
        final String[] arguments = {
                "-" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-10A-01D-1931-08-chr20-chr21.denoisedCR.tsv",
                "-" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-10A-01D-1931-08-chr20-chr21.allelicCounts.tsv",
                "-" + CopyNumberStandardArgument.OUTPUT_PREFIX_SHORT_NAME, "TCGA-05-4389-10A-01D-1931-08-chr20-chr21",
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "/home/slee/working/gatk",
                "-" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testWGS() {
        final String[] arguments = {
                "-" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-01A-01D-1931-08.denoisedCR.tsv",
                "-" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-01A-01D-1931-08.allelicCounts.tsv",
                "-" + CopyNumberStandardArgument.OUTPUT_PREFIX_SHORT_NAME, "TCGA-05-4389-01A-01D-1931-08",
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "/home/slee/working/gatk",
                "-" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testWGSNormal() {
        final String[] arguments = {
//                "-" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-10A-01D-1931-08.denoisedCR.tsv",
                "-" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-10A-01D-1931-08.allelicCounts.tsv",
                "-" + CopyNumberStandardArgument.OUTPUT_PREFIX_SHORT_NAME, "TCGA-05-4389-10A-01D-1931-08",
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "/home/slee/working/gatk",
                "-" + ModelSegments.MAX_NUM_SMOOTHING_ITERATIONS_LONG_NAME, "0",
                "-" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testWESHCC() {
        final String[] arguments = {
                "-" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, "/home/slee/working/gatk/hcc1143_T_clean.denoisedCR.tsv",
                "-" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_SHORT_NAME, "/home/slee/working/gatk/hcc1143_T_clean.allelicCounts.tsv",
                "-" + CopyNumberStandardArgument.OUTPUT_PREFIX_SHORT_NAME, "hcc1143_T_clean",
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "/home/slee/working/gatk",
                "-" + ModelSegments.GENOTYPING_P_VALUE_THRESHOLD_SHORT_NAME, "1",
                "-" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
    }
}