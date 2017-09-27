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
    public void testWGS() {
        final String[] arguments = {
                "-" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-01A-01D-1931-08.chr20-chr21.denoisedCR.tsv",
                "-" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-01A-01D-1931-08.chr20-chr21.allelicCounts.tsv",
                "-" + CopyNumberStandardArgument.OUTPUT_PREFIX_SHORT_NAME, "TCGA-05-4389-01A-01D-1931-08.chr20-chr21",
                "-" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
        };
        runCommandLine(arguments);
    }

//    @Test
//    public void testWES() {
//        final String[] arguments = {
//                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-01A-01D-1265-08-gc-corrected.tn.tsv",
//                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-01A-01D-1265-08-gc-corrected.seg",
//                "-" + StandardArgumentDefinitions.VERBOSITY_NAME, "DEBUG"
//        };
//        runCommandLine(arguments);
//    }
//
//    @Test
//    public void testWGS() {
//        final String[] arguments = {
//                "-" + StandardArgumentDefinitions.INPUT_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-10A-01D-1931-08.coverage.tsv.raw_cov.hdf5.tn.tsv",
//                "-" + StandardArgumentDefinitions.OUTPUT_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-10A-01D-1931-08.coverage.tsv.raw_cov.hdf5.seg",
//                "-" + StandardArgumentDefinitions.VERBOSITY_NAME, "DEBUG"
//        };
//        runCommandLine(arguments);
//    }
//
//    @Test
//    public void testWESAllelic() {
//        final String[] arguments = {
//                "-" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-01A-01D-1265-08-gc-corrected.tn.tsv",
//                "-" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-01A-01D-1265-08.ac.tsv",
//                "-" + CopyNumberStandardArgument.OUTPUT_PREFIX_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-01A-01D-1265-08",
//                "-" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
//        };
//        runCommandLine(arguments);
//    }
//
//    @Test
//    public void testWESAllelicNormal() {
//        final String[] arguments = {
//                "-" + CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-10A-01D-1265-08-gc-corrected.tn.tsv",
//                "-" + CopyNumberStandardArgument.ALLELIC_COUNTS_FILE_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-10A-01D-1265-08.ac.tsv",
//                "-" + CopyNumberStandardArgument.OUTPUT_PREFIX_SHORT_NAME, "/home/slee/working/gatk/TCGA-05-4389-10A-01D-1265-08",
//                "-" + StandardArgumentDefinitions.VERBOSITY_NAME, "INFO"
//        };
//        runCommandLine(arguments);
//    }
}