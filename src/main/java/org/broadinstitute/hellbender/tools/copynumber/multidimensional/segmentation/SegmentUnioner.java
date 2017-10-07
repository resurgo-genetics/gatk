package org.broadinstitute.hellbender.tools.copynumber.multidimensional.segmentation;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.collections4.ListUtils;
import org.apache.commons.math3.util.FastMath;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.exceptions.GATKException;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.allelic.segmentation.AlleleFractionSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.coverage.copyratio.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.coverage.segmentation.CopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.formats.collections.LocatableCollection;
import org.broadinstitute.hellbender.tools.copynumber.utils.segmentation.KernelSegmenter;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.util.*;
import java.util.function.BiFunction;
import java.util.stream.Collectors;
import java.util.stream.Stream;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
final class SegmentUnioner {
    private static final Logger logger = LogManager.getLogger(SegmentUnioner.class);

    //some hardcoded values to use for merging spurious middles via kernel segmentation
    private static final BiFunction<CopyRatio, CopyRatio, Double> COPY_RATIO_KERNEL = (x, y) ->
            x.getLog2CopyRatioValue() * y.getLog2CopyRatioValue();
    private static final double ALLELE_FRACTION_KERNEL_VARIANCE = 0.01;
    private static final BiFunction<AllelicCount, AllelicCount, Double> ALLELE_FRACTION_KERNEL = (x, y) ->
            FastMath.exp(-(x.getAlternateAlleleFraction() - y.getAlternateAlleleFraction()) *
                    (x.getAlternateAlleleFraction() - y.getAlternateAlleleFraction()) /
                    (2. * ALLELE_FRACTION_KERNEL_VARIANCE));
    private static final int KERNEL_APPROXIMATION_DIMENSION = 20;
    private static final int WINDOW_SIZE_DIVISOR = 10;

    private final OverlapDetector<CopyRatio> copyRatioMidpointOverlapDetector;
    private final OverlapDetector<AllelicCount> allelicCountOverlapDetector;
    private final List<SimpleInterval> unionedSegments;

    SegmentUnioner(final CopyRatioSegmentCollection copyRatioSegments,
                   final CopyRatioCollection denoisedCopyRatios,
                   final AlleleFractionSegmentCollection alleleFractionSegments,
                   final AllelicCountCollection allelicCounts) {
        //for determining membership of copy-ratio intervals in unioned segments,
        //we must map copy-ratio intervals to their midpoints so that each will be
        //uniquely contained in a single segment
        copyRatioMidpointOverlapDetector = denoisedCopyRatios.getMidpointOverlapDetector();
        allelicCountOverlapDetector = allelicCounts.getOverlapDetector();
        unionedSegments = constructUnionedSegments(
                copyRatioSegments, denoisedCopyRatios, copyRatioMidpointOverlapDetector,
                alleleFractionSegments, allelicCountOverlapDetector);
        logger.info(unionedSegments.stream().map(copyRatioMidpointOverlapDetector::getOverlaps).mapToDouble(Set::size).sum());
        logger.info(unionedSegments.stream().map(allelicCountOverlapDetector::getOverlaps).mapToDouble(Set::size).sum());
    }

    List<CRAFSegment> constructUnionedCRAFSegments() {
        return unionedSegments.stream()
                .map(s -> new CRAFSegment(
                        s,
                        copyRatioMidpointOverlapDetector.getOverlaps(s).stream()
                                .map(CopyRatio::getLog2CopyRatioValue)
                                .collect(Collectors.toList()),
                        new ArrayList<>(allelicCountOverlapDetector.getOverlaps(s))))
                .collect(Collectors.toList());
    }
    
    /**
     * Returns segments derived from the union of copy-ratio and allele-fraction segments.
     * First, all breakpoints from both sets of segments are combined to form new segments.
     * Spurious segments (i.e., segments containing only copy-ratio intervals that are created by
     * allele-fraction breakpoints and are not present in the original set of copy-ratio segments, and vice versa)
     * at the starts and ends of the original copy-ratio segments are then remerged to the right and left, respectively;
     * spurious segments introduced within the original copy-ratio or allele-fraction segments are merged
     * with adjacent segments by removing the less favorable breakpoint according to
     * the appropriate {@link KernelSegmenter}.  Finally, the segments are trimmed by {@link #trimSegments}.
     */
    private static List<SimpleInterval> constructUnionedSegments(final CopyRatioSegmentCollection copyRatioSegments,
                                                                 final CopyRatioCollection denoisedCopyRatios,
                                                                 final OverlapDetector<CopyRatio> copyRatioMidpointOverlapDetector,
                                                                 final AlleleFractionSegmentCollection alleleFractionSegments,
                                                                 final OverlapDetector<AllelicCount> allelicCountOverlapDetector) {
        final SortedMap<String, List<Breakpoint>> breakpointsByContig = collectBreakpointsByContig(copyRatioSegments, alleleFractionSegments);
        final List<SimpleInterval> untrimmedSegments = constructUntrimmedSegments(copyRatioMidpointOverlapDetector, allelicCountOverlapDetector, breakpointsByContig);
        logger.info(String.format("%d untrimmed segments created...", untrimmedSegments.size()));
        logger.info(untrimmedSegments.stream().map(copyRatioMidpointOverlapDetector::getOverlaps).mapToDouble(Set::size).sum());
        logger.info(untrimmedSegments.stream().map(allelicCountOverlapDetector::getOverlaps).mapToDouble(Set::size).sum());

        //merge spurious segments containing only copy-ratio intervals that were created by allele-fraction breakpoints
        logger.info("Merging spurious copy-ratio segments using default linear kernel...");
        final List<SimpleInterval> spuriousCopyRatioMergedSegments = mergeSpuriousSegments(
                untrimmedSegments, copyRatioSegments, copyRatioMidpointOverlapDetector, allelicCountOverlapDetector, COPY_RATIO_KERNEL);
        logger.info(String.format("%d segments remain after merging spurious copy-ratio segments...", spuriousCopyRatioMergedSegments.size()));
        logger.info(spuriousCopyRatioMergedSegments.stream().map(copyRatioMidpointOverlapDetector::getOverlaps).mapToDouble(Set::size).sum());
        logger.info(spuriousCopyRatioMergedSegments.stream().map(allelicCountOverlapDetector::getOverlaps).mapToDouble(Set::size).sum());

        //merge spurious segments containing only allelic-count sites that were created by copy-ratio breakpoints
        logger.info(String.format("Merging spurious allele-fraction segments using default Gaussian kernel (variance = %4.2f)...", ALLELE_FRACTION_KERNEL_VARIANCE));
        final List<SimpleInterval> spuriousAlleleFractionMergedSegments = mergeSpuriousSegments(
                spuriousCopyRatioMergedSegments, alleleFractionSegments, allelicCountOverlapDetector, copyRatioMidpointOverlapDetector, ALLELE_FRACTION_KERNEL);
        logger.info(String.format("%d segments remain after merging spurious allele-fraction segments...", spuriousAlleleFractionMergedSegments.size()));
        logger.info(spuriousAlleleFractionMergedSegments.stream().map(copyRatioMidpointOverlapDetector::getOverlaps).mapToDouble(Set::size).sum());
        logger.info(spuriousAlleleFractionMergedSegments.stream().map(allelicCountOverlapDetector::getOverlaps).mapToDouble(Set::size).sum());

        logger.info("Trimming combined segments...");
        //for trimming segments in the final step, we need an overlap detector built from the full copy-ratio intervals
        final OverlapDetector<CopyRatio> copyRatioOverlapDetector = denoisedCopyRatios.getOverlapDetector();
        return spuriousAlleleFractionMergedSegments.stream()
                .map(s -> trimSegments(s, copyRatioOverlapDetector, allelicCountOverlapDetector))
                .filter(Objects::nonNull)
                .collect(Collectors.toList());
    }

    /**
     * Returns a map of contig -> combined breakpoints from copy-ratio and allele-fraction segments on that contig.
     */
    private static SortedMap<String, List<Breakpoint>> collectBreakpointsByContig(final CopyRatioSegmentCollection copyRatioSegments,
                                                                                  final AlleleFractionSegmentCollection alleleFractionSegments) {
        return ListUtils.union(copyRatioSegments.getIntervals(), alleleFractionSegments.getIntervals()).stream()
                .map(s -> Arrays.asList(
                        new Breakpoint(s.getContig(), s.getStart(), BreakpointType.START),
                        new Breakpoint(s.getContig(), s.getEnd(), BreakpointType.END)))
                .flatMap(Collection::stream)
                .sorted((b1, b2) -> LocatableCollection.LEXICOGRAPHICAL_ORDER_COMPARATOR.compare(b1.interval, b2.interval))
                .collect(Collectors.groupingBy(Breakpoint::getContig, TreeMap::new, Collectors.toList()));
    }

    /**
     * Given breakpointsByContig map, constructs a list of untrimmed segments (i.e., segments that are directly
     * adjacent to each other; segment boundaries at each breakpoint are determined by whether that breakpoint
     * was a start or an end for the corresponding original copy-ratio/allele-fraction segment),
     * then returns only those untrimmed segments that are non-empty.
     */
    private static List<SimpleInterval> constructUntrimmedSegments(final OverlapDetector<CopyRatio> copyRatioMidpointOverlapDetector,
                                                                   final OverlapDetector<AllelicCount> allelicCountOverlapDetector,
                                                                   final SortedMap<String, List<Breakpoint>> breakpointsByContig) {
        final List<SimpleInterval> segments = new ArrayList<>();
        for (final String contig : breakpointsByContig.keySet()) {
            final List<Breakpoint> breakpoints = breakpointsByContig.get(contig);
            int start = breakpoints.get(0).getSite();
            for (int i = 1; i < breakpoints.size(); i++) {
                final Breakpoint breakpoint = breakpoints.get(i);

                //adjust segment boundaries if next breakpoint is at the same location
                // or according to breakpoint type (this prevents, e.g., allelic-count sites
                //that originally started allele-fraction segments from being stranded)
                final boolean isNextBreakpointSameLocation = i < breakpoints.size() - 1 &&
                        breakpoints.get(i + 1).getSite() == breakpoint.getSite();
                final int adjustment =
                        breakpoint.type == BreakpointType.START ||
                                (breakpoint.type == BreakpointType.END && isNextBreakpointSameLocation)
                                ? 1 : 0;
                final int end = breakpoint.getSite() - adjustment;

                if (end < start) {  //this could happen if there are adjacent breakpoints
                    continue;
                }

                final SimpleInterval segment = new SimpleInterval(contig, start, end);
                if (copyRatioMidpointOverlapDetector.overlapsAny(segment) || allelicCountOverlapDetector.overlapsAny(segment)) {
                    segments.add(segment);
                }
                start = end + 1;
            }
        }
        return segments;
    }

    /**
     * Given combined segments, returns a list with spurious segments---which contain only points
     * of the first type that were created during breakpoint union by breakpoints arising from
     * points of the second type---remerged.
     */
    private static <FirstTypeSegment extends Locatable, FirstType extends Locatable, SecondType extends Locatable> List<SimpleInterval> mergeSpuriousSegments(
            final List<SimpleInterval> combinedSegments,
            final LocatableCollection<FirstTypeSegment> originalFirstTypeSegments,
            final OverlapDetector<FirstType> firstTypeOverlapDetector,
            final OverlapDetector<SecondType> secondTypeOverlapDetector,
            final BiFunction<FirstType, FirstType, Double> kernelForFirstTypeSegmentation) {
        //get starts and ends of original first-type segments
        final Set<SimpleInterval> originalFirstTypeSegmentStarts =
                originalFirstTypeSegments.getIntervals().stream()
                        .map(s -> new SimpleInterval(s.getContig(), s.getStart(), s.getStart()))
                        .collect(Collectors.toSet());
        final Set<SimpleInterval> originalFirstTypeSegmentEnds =
                originalFirstTypeSegments.getIntervals().stream()
                        .map(s -> new SimpleInterval(s.getContig(), s.getEnd(), s.getEnd()))
                        .collect(Collectors.toSet());

        final List<SimpleInterval> mergedSegments = new ArrayList<>();
        final ListIterator<SimpleInterval> segmentsIter = combinedSegments.listIterator();
        while (segmentsIter.hasNext()) {
            final SimpleInterval segment = segmentsIter.next();
            //do not remerge segments containing points of the other type or segments from the original segmentation
            if (secondTypeOverlapDetector.overlapsAny(segment) || originalFirstTypeSegments.getIntervals().contains(segment)) {
                mergedSegments.add(segment);
                continue;
            }
            final SimpleInterval segmentStart =
                    new SimpleInterval(segment.getContig(), segment.getStart(), segment.getStart());
            final SimpleInterval segmentEnd =
                    new SimpleInterval(segment.getContig(), segment.getEnd(), segment.getEnd());
            if (originalFirstTypeSegmentStarts.contains(segmentStart) && !originalFirstTypeSegmentEnds.contains(segmentEnd)) {
                //remerge segments introduced at starts to the right
                final SimpleInterval nextSegment = segmentsIter.next();
                mergedSegments.add(mergeSegments(segment, nextSegment));
            } else if (!originalFirstTypeSegmentStarts.contains(segmentStart) && originalFirstTypeSegmentEnds.contains(segmentEnd)) {
                //remerge segments introduced at ends to the left
                final int previousIndex = mergedSegments.size() - 1;
                final SimpleInterval previousSegment = mergedSegments.get(previousIndex);
                mergedSegments.set(previousIndex, mergeSegments(previousSegment, segment));
            } else if (!originalFirstTypeSegmentStarts.contains(segmentStart) && !originalFirstTypeSegmentEnds.contains(segmentEnd)) {
                //remerge segments introduced in the middle by examining adjacent segments
                logger.debug(String.format("Merging spurious middle segment: %s (number of data points = %d)",
                        segment, firstTypeOverlapDetector.getOverlaps(segment).size()));
                final int previousIndex = mergedSegments.size() - 1;
                final SimpleInterval previousSegment = mergedSegments.get(previousIndex);
                final SimpleInterval nextSegment = combinedSegments.get(combinedSegments.indexOf(segment) + 1);

                //if either adjacent segment does not contain points of the first type, then do not merge with adjacent segments
                if (firstTypeOverlapDetector.getOverlaps(previousSegment).isEmpty() || firstTypeOverlapDetector.getOverlaps(nextSegment).isEmpty()) {
                    logger.debug("An adjacent segment was empty, not merging spurious middle segment...");
                    mergedSegments.add(segment);
                    continue;
                }

                //get points of the first type from current and adjacent segments
                final List<FirstType> firstTypePoints = Stream.of(previousSegment, segment, nextSegment)
                        .flatMap(s -> firstTypeOverlapDetector.getOverlaps(s).stream())
                        .sorted(LocatableCollection.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                        .collect(Collectors.toList());

                //create kernel segmenter for points of first type contained in all three segments (left + middle + right)
                final KernelSegmenter<FirstType> segmenter = new KernelSegmenter<>(firstTypePoints);

                //test for single changepoint with one window size and no penalty
                final int bestChangepointIndex = segmenter.findChangepoints(
                        1,
                        kernelForFirstTypeSegmentation,
                        Math.min(KERNEL_APPROXIMATION_DIMENSION, firstTypePoints.size()),
                        Collections.singletonList(Math.max(firstTypePoints.size() / WINDOW_SIZE_DIVISOR, 1)),    //try to pick a relatively big window size, but we are only guaranteed to have at least 3 points)
                        0., 0.,                     //no penalty guarantees one changepoint
                        KernelSegmenter.ChangepointSortOrder.BACKWARD_SELECTION).get(0);
                final FirstType bestChangepoint = firstTypePoints.get(bestChangepointIndex);
                logger.debug(String.format("Best changepoint: %s", bestChangepoint));

                //if changepoint falls within middle (i.e., spurious) segment,
                //then create new breakpoint there and merge resulting segments to the left and right
                //else if changepoint falls in left (right) segment, then merge middle segment to the right (left)
                if (firstTypeOverlapDetector.getOverlaps(segment).contains(bestChangepoint)) {
                    final FirstType pointAfterBestChangepoint = firstTypePoints.get(bestChangepointIndex + 1);
                    final SimpleInterval newLeftSegment = new SimpleInterval(segment.getContig(), previousSegment.getStart(), bestChangepoint.getEnd());
                    final SimpleInterval newRightSegment = new SimpleInterval(segment.getContig(), pointAfterBestChangepoint.getStart(), nextSegment.getEnd());
                    logger.debug(String.format("Best changepoints found in current segment, creating new segments: %s, %s", newLeftSegment, newRightSegment));
                    mergedSegments.set(previousIndex, newLeftSegment);
                    mergedSegments.add(newRightSegment);
                    segmentsIter.next();
                } else if (firstTypeOverlapDetector.getOverlaps(previousSegment).contains(bestChangepoint)) {
                    logger.debug("Best changepoint found in adjacent segment to left, merging middle segment to the right...");
                    mergedSegments.add(mergeSegments(segment, nextSegment));
                    segmentsIter.next();
                } else {
                    logger.debug("Best changepoint found in adjacent segment to right, merging middle segment to the left...");
                    mergedSegments.set(previousIndex, mergeSegments(previousSegment, segment));
                }
            }
        }
        return mergedSegments;
    }

    /**
     * Given a segment and collections of copy ratios and allelic counts, returns a trimmed segment produced by
     * removing the empty portions at the start and the end of the original segment that do not overlap the
     * copy ratios and allelic counts that overlap with the original segment.
     * If this procedure would remove the entire segment, {@code null} is returned instead.
     */
    private static SimpleInterval trimSegments(final SimpleInterval segment,
                                               final OverlapDetector<CopyRatio> copyRatioOverlapDetector,
                                               final OverlapDetector<AllelicCount> allelicCountOverlapDetector) {
        final Set<CopyRatio> copyRatiosOverlappingSegment = copyRatioOverlapDetector.getOverlaps(segment);
        final Set<AllelicCount> allelicCountsInSegment = allelicCountOverlapDetector.getOverlaps(segment);

        final int numCopyRatiosOverlappingSegment = copyRatiosOverlappingSegment.size();
        final int numAllelicCountsInSegment = allelicCountsInSegment.size();

        if (numCopyRatiosOverlappingSegment == 0 && numAllelicCountsInSegment == 0) {
            return null;
        }

        int start = segment.getStart();
        int end = segment.getEnd();

        if (numCopyRatiosOverlappingSegment == 0 && numAllelicCountsInSegment > 0) {
            //if there are no copy-ratio intervals overlapping segment, use allelic counts to determine trimmed segment
            start = allelicCountsInSegment.stream().mapToInt(AllelicCount::getStart).min().getAsInt();
            end = allelicCountsInSegment.stream().mapToInt(AllelicCount::getEnd).max().getAsInt();
        } else if (numCopyRatiosOverlappingSegment > 0) {
            //if segment start does not fall within first copy-ratio interval, use start of first copy-ratio interval as start of trimmed segment
            start = Math.max(start, copyRatiosOverlappingSegment.stream().mapToInt(CopyRatio::getStart).min().getAsInt());
            //if segment end does not fall within last copy-ratio interval, use end of last copy-ratio interval as end of trimmed segment
            end = Math.min(end, copyRatiosOverlappingSegment.stream().mapToInt(CopyRatio::getEnd).max().getAsInt());
            //if there are also allelic counts within segment, check to see if they give a larger trimmed segment
            if (numAllelicCountsInSegment > 0) {
                start = Math.min(start, allelicCountsInSegment.stream().mapToInt(AllelicCount::getStart).min().getAsInt());
                end = Math.max(end, allelicCountsInSegment.stream().mapToInt(AllelicCount::getEnd).max().getAsInt());
            }
        }
        if (start < segment.getStart() || end > segment.getEnd() || end < start) {
            throw new GATKException.ShouldNeverReachHereException("Something went wrong in trimming interval.");
        }
        return new SimpleInterval(segment.getContig(), start, end);
    }

    private static SimpleInterval mergeSegments(final SimpleInterval segment1,
                                                final SimpleInterval segment2) {
        Utils.validateArg(segment1.getContig().equals(segment2.getContig()),
                String.format("Cannot join segments %s and %s on different chromosomes.", segment1.toString(), segment2.toString()));
        final int start = Math.min(segment1.getStart(), segment2.getStart());
        final int end = Math.max(segment1.getEnd(), segment2.getEnd());
        return new SimpleInterval(segment1.getContig(), start, end);
    }

    private enum BreakpointType {
        START, END
    }

    private static final class Breakpoint {
        private final SimpleInterval interval;
        private final BreakpointType type;

        Breakpoint(final String contig,
                   final int site,
                   final BreakpointType type) {
            interval = new SimpleInterval(contig, site, site);
            this.type = type;
        }

        BreakpointType getType() {
            return type;
        }

        String getContig() {
            return interval.getContig();
        }

        int getSite() {
            return interval.getStart();
        }
    }
}