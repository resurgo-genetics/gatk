package org.broadinstitute.hellbender.tools.exome;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.stream.Collectors;

/**
 * Represents copy-number data from a single individual for exome analysis.  Contains coverage data from targets
 * and ref/alt allele counts at normal germline het SNP sites.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class Genome {
    private final TargetCollection<ReadCountRecord.SingleSampleRecord> targets;
    private final TargetCollection<AllelicCount> snps;
    private final String sampleName;

    /**
     * Constructs a genome from lists containing log_2 target-coverage and SNP-allele-count data.
     * @param targets       list of log_2 target coverages, cannot be {@code null}
     * @param snps          list of SNP allele counts, cannot be {@code null}
     */
    public <T extends ReadCountRecord> Genome(final List<T> targets, final List<AllelicCount> snps, final String sampleName) {
        Utils.nonNull(targets);
        Utils.nonNull(snps);
        Utils.nonNull(sampleName);
        this.targets = new HashedListTargetCollection<>(targets.stream().map(ReadCountRecord::asSingleSampleRecord).collect(Collectors.toList()));
        this.snps = new HashedListTargetCollection<>(snps);
        this.sampleName = sampleName;
    }

    public Genome(final CopyRatioCollection denoisedCopyRatios, final org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollection allelicCounts) {
        Utils.nonNull(denoisedCopyRatios);
        Utils.nonNull(allelicCounts);
        Utils.validateArg(denoisedCopyRatios.getSampleName().equals(allelicCounts.getSampleName()),
                "Sample names from copy ratios and allelic counts must match.");
        this.targets = new HashedListTargetCollection<>(denoisedCopyRatios.getRecords().stream()
                .map(cr -> new ReadCountRecord.SingleSampleRecord(new Target(cr.getInterval()), cr.getLog2CopyRatioValue()))
                .collect(Collectors.toList()));
        this.snps = new HashedListTargetCollection<>(allelicCounts.getRecords().stream()
                .map(ac -> new AllelicCount(ac.getInterval(), ac.getRefReadCount(), ac.getAltReadCount()))
                .collect(Collectors.toList()));
        sampleName = denoisedCopyRatios.getSampleName();
    }

    /**
     * Constructs a genome from files containing log_2 target-coverage and SNP-allele-count data.
     * @param tangentNormalizedCoverageFile     log_2 target-coverage file
     * @param snpFile                           SNP-allele-count file
     */
    public Genome(final File tangentNormalizedCoverageFile, final File snpFile) {
        sampleName = ReadCountCollectionUtils.getSampleNameForCLIsFromReadCountsFile(tangentNormalizedCoverageFile);
        try {
            targets = new HashedListTargetCollection<>(ReadCountCollectionUtils.parse(tangentNormalizedCoverageFile).records()
                    .stream().map(ReadCountRecord::asSingleSampleRecord).collect(Collectors.toList()));
            snps = new HashedListTargetCollection<>(new AllelicCountCollection(snpFile).getCounts());
        } catch (final IOException e) {
            throw new UserException.BadInput("Could not read normalized coverage file");
        }
    }

    public final TargetCollection<ReadCountRecord.SingleSampleRecord> getTargets() {  return targets; }

    public final TargetCollection<AllelicCount> getSNPs() {  return snps; }

    public final String getSampleName() {   return sampleName;  }
}
