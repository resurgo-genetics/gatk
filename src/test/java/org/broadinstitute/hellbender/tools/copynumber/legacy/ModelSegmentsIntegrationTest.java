package org.broadinstitute.hellbender.tools.copynumber.legacy;

import org.broadinstitute.hellbender.CommandLineProgramTest;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.CopyNumberStandardArgument;
import org.junit.Test;

/**
 * Created by slee on 9/6/17.
 */
public class ModelSegmentsIntegrationTest extends CommandLineProgramTest {
    @Test
    public void testWGSChr20Chr21() {
        final String[] arguments = {
                "-" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-01A-01D-1931-08-chr20-chr21.denoisedCR.tsv",
                "-" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-01A-01D-1931-08-chr20-chr21.allelicCounts.tsv",
                "-" + CopyNumberStandardArgument.OUTPUT_PREFIX_SHORT_NAME, "TCGA-05-4389-01A-01D-1931-08-chr20-chr21",
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "/home/slee/working/gatk",
                "--" + ModelSegments.NUM_SAMPLES_COPY_RATIO_LONG_NAME, "50",
                "--" + ModelSegments.NUM_BURN_IN_COPY_RATIO_LONG_NAME, "10",
                "--" + ModelSegments.NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME, "50",
                "--" + ModelSegments.NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME, "10",
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
                "--" + ModelSegments.NUM_SAMPLES_COPY_RATIO_LONG_NAME, "50",
                "--" + ModelSegments.NUM_BURN_IN_COPY_RATIO_LONG_NAME, "10",
                "--" + ModelSegments.NUM_SAMPLES_ALLELE_FRACTION_LONG_NAME, "50",
                "--" + ModelSegments.NUM_BURN_IN_ALLELE_FRACTION_LONG_NAME, "10",
                "-" + StandardArgumentDefinitions.VERBOSITY_NAME, "DEBUG"
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
                "-" + StandardArgumentDefinitions.VERBOSITY_NAME, "DEBUG"
        };
        runCommandLine(arguments);
    }

    @Test
    public void testWGSNormal() {
        final String[] arguments = {
                "-" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-10A-01D-1931-08.denoisedCR.tsv",
                "-" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-10A-01D-1931-08.allelicCounts.tsv",
                "-" + CopyNumberStandardArgument.OUTPUT_PREFIX_SHORT_NAME, "TCGA-05-4389-10A-01D-1931-08",
                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "/home/slee/working/gatk",
                "-" + ModelSegments.MAXIMUM_NUMBER_OF_SEGMENTS_PER_CHROMOSOME_LONG_NAME, "1",
                "-" + StandardArgumentDefinitions.VERBOSITY_NAME, "DEBUG"
        };
        runCommandLine(arguments);
    }
}