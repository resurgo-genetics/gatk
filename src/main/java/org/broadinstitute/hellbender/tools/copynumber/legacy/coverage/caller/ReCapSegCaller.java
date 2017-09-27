package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.caller;

import org.apache.commons.math3.stat.descriptive.moment.Mean;
import org.apache.commons.math3.stat.descriptive.moment.StandardDeviation;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioSegment;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioSegmentCollection;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.stream.Collectors;

/**
 * <p>This caller mimics the legacy ReCapSeg Caller that was originally implemented in ReCapSeg v1.4.5.0.</p>
 *
 * <p>There is a small difference.  The python code was using the same algorithm as intersectBed, which was causing it to drop
 *  the first target of each segment in calculations of the copy neutral targets.  The code here (targets.targets(s)) does
 *  not do this.  This difference in the two codebases can cause a slight difference in the T calculation.  Hence, the
 *  results of this code and the python code will not be exactly the same, but will be
 *  very close.  A fix (to make this code match the python) has been deemed unworthy of our time.</p>
 *
 *  Note: Assumes that the input @{link ReadCountCollection} contains only one sample -- calls will be made for only one
 *        counts column and other columns ignored.
 *
 */
public final class ReCapSegCaller {
    private static final Logger logger = LogManager.getLogger(ReCapSegCaller.class);

    //bounds on log_2 coverage for high-confidence neutral segments
    private static final double COPY_NEUTRAL_CUTOFF = 0.1;
    // Number of standard deviations before assuming that a target was an outlier in a segment
    private static final double Z_THRESHOLD = 2;

    private final CopyRatioCollection denoisedCopyRatios;
    private final CopyRatioSegmentCollection copyRatioSegments;
    private final LinkedHashMap<CopyRatioSegment, List<CopyRatio>> segmentToCopyRatiosMap;

    /**
     * @param denoisedCopyRatios in log2 space
     */
    public ReCapSegCaller(final CopyRatioCollection denoisedCopyRatios,
                          final CopyRatioSegmentCollection copyRatioSegments) {
        this.denoisedCopyRatios = Utils.nonNull(denoisedCopyRatios);
        this.copyRatioSegments = Utils.nonNull(copyRatioSegments);
        Utils.validateArg(denoisedCopyRatios.getSampleName().equals(copyRatioSegments.getSampleName()),
                "Denoised copy ratios and copy-ratio segments do not have the same sample name.");
        segmentToCopyRatiosMap = constructSegmentToCopyRatiosMap(denoisedCopyRatios, copyRatioSegments);
    }

    private static LinkedHashMap<CopyRatioSegment, List<CopyRatio>> constructSegmentToCopyRatiosMap(final CopyRatioCollection denoisedCopyRatios,
                                                                                                    final CopyRatioSegmentCollection copyRatioSegments) {
        final LinkedHashMap<CopyRatioSegment, List<CopyRatio>> segmentToCopyRatiosMap = new LinkedHashMap<>();
        int index = 0;
        for (final CopyRatioSegment segment : copyRatioSegments.getRecords()) {
            final int numPoints = segment.getNumPoints();
            final int firstPointStart = denoisedCopyRatios.getRecords().get(index).getStart();
            final int lastPointEnd = denoisedCopyRatios.getRecords().get(index + numPoints - 1).getEnd();
            if (!(firstPointStart == segment.getStart() && lastPointEnd == segment.getEnd())) {
                throw new IllegalArgumentException("Denoised copy ratios and copy-ratio segments are not consistent.");
            }
            segmentToCopyRatiosMap.put(segment, denoisedCopyRatios.getRecords().subList(index, index + numPoints));
            index += numPoints;
        }
        if (index != denoisedCopyRatios.getRecords().size()) {
            throw new IllegalArgumentException("Denoised copy ratios and copy-ratio segments are not consistent.");
        }
        return segmentToCopyRatiosMap;
    }

    private double calculateT() {
        //Get the segments that are likely copy neutral.
        // Math.abs removed to mimic python...
        final List<CopyRatioSegment> copyNeutralSegments = segmentToCopyRatiosMap.keySet().stream()
                .filter(s -> s.getMeanLog2CopyRatio() < COPY_NEUTRAL_CUTOFF).collect(Collectors.toList());

        // Get the targets that correspond to the copyNeutralSegments... note that individual targets, due to noise,
        //  can be far away from copy neutral
        final double[] copyNeutralTargetsCopyRatio = copyNeutralSegments.stream()
                .flatMap(s -> segmentToCopyRatiosMap.get(s).stream())
                .mapToDouble(CopyRatio::getLog2CopyRatioValue).toArray();

        final double meanCopyNeutralTargets = new Mean().evaluate(copyNeutralTargetsCopyRatio);
        final double sigmaCopyNeutralTargets = new StandardDeviation().evaluate(copyNeutralTargetsCopyRatio);

        // Now we filter outliers by only including those w/in 2 standard deviations.
        final double [] filteredCopyNeutralTargetsCopyRatio = Arrays.stream(copyNeutralTargetsCopyRatio)
                .filter(c -> Math.abs(c - meanCopyNeutralTargets) < sigmaCopyNeutralTargets * Z_THRESHOLD).toArray();

        return new StandardDeviation().evaluate(filteredCopyNeutralTargetsCopyRatio);
    }

    public CalledCopyRatioSegmentCollection makeCalls() {
        final double t = calculateT();

        logger.info("Running caller that mimics the ReCapSeg 1.4.5.0 (python) caller.");
        // Log some information about thresholds chosen for the segments.
        logger.info(String.format("Copy neutral (log2CR space) [%.4f, %.4f]", -t, t));
        logger.info(String.format("Copy neutral (CR space) [%.4f, %.4f]", Math.pow(2, -t), Math.pow(2, t)));

        final Set<CopyRatioSegment> segments = segmentToCopyRatiosMap.keySet();
        final List<CalledCopyRatioSegment> calledSegments = new ArrayList<>(segments.size());
        for (final CopyRatioSegment segment : segments) {
            if (segment.getMeanLog2CopyRatio() < -t) {
                calledSegments.add(new CalledCopyRatioSegment(segment, CalledCopyRatioSegment.Call.DELETION_CALL));
            } else if (segment.getMeanLog2CopyRatio() > t) {
                calledSegments.add(new CalledCopyRatioSegment(segment, CalledCopyRatioSegment.Call.AMPLIFICATION_CALL));
            } else {
                calledSegments.add(new CalledCopyRatioSegment(segment, CalledCopyRatioSegment.Call.NEUTRAL_CALL));
            }
        }

        return new CalledCopyRatioSegmentCollection(copyRatioSegments.getSampleName(), calledSegments);
    }
}
