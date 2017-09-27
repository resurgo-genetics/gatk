package org.broadinstitute.hellbender.tools.copynumber.legacy.plotting;

import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.SAMSequenceRecord;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.barclay.help.DocumentedFeature;
import org.broadinstitute.hellbender.cmdline.CommandLineProgram;
import org.broadinstitute.hellbender.cmdline.StandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.CopyNumberStandardArgument;
import org.broadinstitute.hellbender.utils.R.RScriptExecutor;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.io.Resource;
import org.broadinstitute.hellbender.utils.reference.ReferenceUtils;

import java.io.File;
import java.util.*;
import java.util.stream.Collectors;

/**
 * Plots segmented coverage results.
 *
 * <p>The order and representation of contigs in plots follows the contig ordering within the required reference sequence dictionary. </p>
 *
 * <h3>Examples</h3>
 *
 * <p>The --output parameter specifies a pre-existing directory.</p>
 *
 * <pre>
 * gatk-launch --javaOptions "-Xmx4g" PlotSegmentedCopyRatio \
 *   --standardizedCopyRatios tumor.standardizedCR.tsv \
 *   --standardizedCopyRatios tumor.denoisedCR.tsv \
 *   -SD ref_fasta.dict \
 *   --output output_file
 * </pre>
 */
@CommandLineProgramProperties(
        summary = "Create plots of denoised copy ratio.",
        oneLineSummary = "Create plots of denoised copy ratio.",
        programGroup = CopyNumberProgramGroup.class
)
@DocumentedFeature
public final class PlotDenoisedCopyRatios extends CommandLineProgram {
    private static final String CNV_PLOTTING_R_LIBRARY = "CNVPlottingLibrary.R";
    private static final String PLOT_DENOISED_COPY_RATIOS_R_SCRIPT = "PlotDenoisedCopyRatios.R";

    private static final String CONTIG_DELIMITER = "CONTIG_DELIMITER";  //used to delimit contig names and lengths passed to the R script
    private static final int DEFAULT_MINIMUM_CONTIG_LENGTH = 1000000;   //can be used to filter out mitochondrial contigs, unlocalized contigs, etc.

    private static final String MINIMUM_CONTIG_LENGTH_LONG_NAME = "minimumContigLength";
    private static final String MINIMUM_CONTIG_LENGTH_SHORT_NAME = "minContigLength";

    @Argument(
            doc = "Input file containing standardized copy-ratio profile (output of DenoiseReadCounts).",
            fullName = CopyNumberStandardArgument.STANDARDIZED_COPY_RATIOS_FILE_FULL_NAME,
            shortName = CopyNumberStandardArgument.STANDARDIZED_COPY_RATIOS_FILE_SHORT_NAME,
            optional = true
    )
    private File inputStandardizedCopyRatiosFile;

    @Argument(
            doc = "Input file containing denoised copy-ratio profile (output of DenoiseReadCounts).",
            fullName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_FULL_NAME,
            shortName = CopyNumberStandardArgument.DENOISED_COPY_RATIOS_FILE_SHORT_NAME,
            optional = true
    )
    private File inputDenoisedCopyRatiosFile;

    @Argument(
            doc = "File containing the reference sequence dictionary (used to determine relative contig lengths). " +
                    "Contigs will be plotted in the order given. " +
                    "Contig names should not include the string \"" + CONTIG_DELIMITER + "\". " +
                    "The tool only considers contigs in the given dictionary for plotting, and " +
                    "data for contigs absent in the dictionary generate only a warning. In other words, you may " +
                    "modify a reference dictionary for use with this tool to include only contigs for which plotting is desired, " +
                    "and sort the contigs to the order in which the plots should display the contigs.",
            shortName = StandardArgumentDefinitions.SEQUENCE_DICTIONARY_SHORT_NAME
    )
    private File sequenceDictionaryFile;

    @Argument(
            doc = "Threshold length (in bp) for contigs to be plotted. " +
                    "Contigs with lengths less than this threshold will not be plotted. " +
                    "This can be used to filter out mitochondrial contigs, unlocalized contigs, etc.",
            fullName =  MINIMUM_CONTIG_LENGTH_LONG_NAME,
            shortName = MINIMUM_CONTIG_LENGTH_SHORT_NAME,
            optional = true
    )
    private int minContigLength = DEFAULT_MINIMUM_CONTIG_LENGTH;

    @Argument(
            doc = "Prefix for output filenames.",
            fullName =  CopyNumberStandardArgument.OUTPUT_PREFIX_LONG_NAME,
            shortName = CopyNumberStandardArgument.OUTPUT_PREFIX_SHORT_NAME
    )
    private String outputPrefix;

    @Argument(
            doc = "Output directory.",
            fullName =  StandardArgumentDefinitions.OUTPUT_LONG_NAME,
            shortName = StandardArgumentDefinitions.OUTPUT_SHORT_NAME
    )
    private String outputDir;

    @Override
    protected Object doWork() {
        checkRegularReadableUserFiles();

        //read input files
        final CopyRatioCollection standardizedCopyRatios = new CopyRatioCollection(inputStandardizedCopyRatiosFile);
        final CopyRatioCollection denoisedCopyRatios = new CopyRatioCollection(inputDenoisedCopyRatiosFile);
        Utils.validateArg(standardizedCopyRatios.getIntervals().equals(denoisedCopyRatios.getIntervals()),
                "Intervals in input files must be identical.");

        //get sample name from input files (consistency check is performed)
        final String sampleName = getSampleName(standardizedCopyRatios, denoisedCopyRatios);

        //load contig names and lengths from the sequence dictionary into a LinkedHashMap
        final SAMSequenceDictionary sequenceDictionary = ReferenceUtils.loadFastaDictionary(sequenceDictionaryFile);
        Utils.validateArg(sequenceDictionary.getSequences().stream().map(SAMSequenceRecord::getSequenceName).noneMatch(n -> n.contains(CONTIG_DELIMITER)),
                String.format("Contig names cannot contain \"%s\".", CONTIG_DELIMITER));
        final Map<String, Integer> contigLengthMap = sequenceDictionary.getSequences().stream()
                .filter(s -> s.getSequenceLength() >= minContigLength)
                .collect(Collectors.toMap(SAMSequenceRecord::getSequenceName, SAMSequenceRecord::getSequenceLength,
                        (c, l) -> {
                            throw new IllegalArgumentException(String.format("Duplicate contig in sequence dictionary: %s", c));
                        },
                        LinkedHashMap::new));
        Utils.validateArg(contigLengthMap.size() > 0,
                "There must be at least one contig above the threshold length in the sequence dictionary.");
        logger.info("Contigs above length threshold: " + contigLengthMap.toString());

        //check that contigs in input files are present in sequence dictionary and that data points are valid given lengths
        validateContigs(inputStandardizedCopyRatiosFile, standardizedCopyRatios, contigLengthMap);
        validateContigs(inputDenoisedCopyRatiosFile, denoisedCopyRatios, contigLengthMap);

        //generate the plots
        final List<String> contigNames = new ArrayList<>(contigLengthMap.keySet());
        final List<Integer> contigLengths = new ArrayList<>(contigLengthMap.values());
        writeDenoisingPlots(sampleName, contigNames, contigLengths);

        return "SUCCESS";
    }

    private void checkRegularReadableUserFiles() {
        IOUtils.canReadFile(inputStandardizedCopyRatiosFile);
        IOUtils.canReadFile(inputDenoisedCopyRatiosFile);
        IOUtils.canReadFile(sequenceDictionaryFile);
        if (!new File(outputDir).exists()) {
            throw new UserException(String.format("Output directory %s does not exist.", outputDir));
        }
    }

    private String getSampleName(final CopyRatioCollection standardizedCopyRatios,
                                 final CopyRatioCollection denoisedCopyRatios) {
        final String standardizedSampleName = standardizedCopyRatios.getSampleName();
        final String denoisedSampleName = denoisedCopyRatios.getSampleName();
        Utils.validateArg(standardizedSampleName.equals(denoisedSampleName),"Sample names in input files must be identical.");
        return standardizedSampleName;
    }

    //validate contig names and lengths
    private void validateContigs(final File copyRatiosFile,
                                 final CopyRatioCollection copyRatios,
                                 final Map<String, Integer> contigLengthMap) {
        final Set<String> contigNamesFromMap = contigLengthMap.keySet();
        final Set<String> contigNames = copyRatios.getRecords().stream().map(CopyRatio::getContig).collect(Collectors.toSet());
        if (!contigNamesFromMap.containsAll(contigNames)) {
            logger.warn(String.format("Contigs present in %s are missing from the sequence dictionary and will not be plotted.", copyRatiosFile.getAbsolutePath()));
        }
        final Map<String, Integer> contigMaxPositionMap = copyRatios.getRecords().stream()
                .filter(i -> contigNamesFromMap.contains(i.getContig()))
                .collect(Collectors.toMap(CopyRatio::getContig, CopyRatio::getEnd, Integer::max));
        contigMaxPositionMap.keySet().forEach(c -> Utils.validateArg(contigMaxPositionMap.get(c) <= contigLengthMap.get(c),
                String.format("Position present in %s exceeds contig length in the sequence dictionary.", copyRatiosFile.getAbsolutePath())));
    }

    /**
     * @param sampleName Sample name derived from input files
     * @param contigNames List containing contig names
     * @param contigLengths List containing contig lengths (same order as contigNames)
     */
    private void writeDenoisingPlots(final String sampleName,
                                     final List<String> contigNames,
                                     final List<Integer> contigLengths) {
        final String contigNamesArg = contigNames.stream().collect(Collectors.joining(CONTIG_DELIMITER));                            //names separated by delimiter
        final String contigLengthsArg = contigLengths.stream().map(Object::toString).collect(Collectors.joining(CONTIG_DELIMITER));  //names separated by delimiter
        final String outputDirArg = addTrailingSlashIfNecessary(outputDir);

        final RScriptExecutor executor = new RScriptExecutor();

        //this runs the R statement "source("CNVPlottingLibrary.R")" before the main script runs
        executor.addScript(new Resource(CNV_PLOTTING_R_LIBRARY, PlotDenoisedCopyRatios.class));
        executor.addScript(new Resource(PLOT_DENOISED_COPY_RATIOS_R_SCRIPT, PlotDenoisedCopyRatios.class));
        //--args is needed for Rscript to recognize other arguments properly
        executor.addArgs("--args",
                "--sample_name=" + sampleName,
                "--standardized_file=" + inputStandardizedCopyRatiosFile,
                "--denoised_file=" + inputDenoisedCopyRatiosFile,
                "--contig_names=" + contigNamesArg,
                "--contig_lengths=" + contigLengthsArg,
                "--output_dir=" + outputDirArg,
                "--output_prefix=" + outputPrefix);
        executor.exec();
    }

    private static String addTrailingSlashIfNecessary(final String outputDir) {
        return outputDir.endsWith(File.separator) ? outputDir : outputDir + File.separator;
    }
}
