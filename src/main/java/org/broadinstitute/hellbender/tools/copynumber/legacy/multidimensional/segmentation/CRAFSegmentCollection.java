package org.broadinstitute.hellbender.tools.copynumber.legacy.multidimensional.segmentation;

import htsjdk.samtools.util.Interval;
import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.apache.commons.collections4.ListUtils;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.segmentation.AlleleFractionSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatio;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.copyratio.CopyRatioCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioKernelSegmenter;
import org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation.CopyRatioSegmentCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.TSVLocatableCollection;
import org.broadinstitute.hellbender.tools.exome.Genome;
import org.broadinstitute.hellbender.tools.exome.SegmentUtils;
import org.broadinstitute.hellbender.tools.exome.SegmentedGenome;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.*;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class CRAFSegmentCollection extends TSVLocatableCollection<CRAFSegment> {
    private static final Logger logger = LogManager.getLogger(CRAFSegmentCollection.class);

    enum CRAFSegmentTableColumn {
        CONTIG,
        START,
        END,
        NUM_POINTS_COPY_RATIO,
        NUM_POINTS_ALLELE_FRACTION,
        MEAN_LOG2_COPY_RATIO;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private static final Function<DataLine, CRAFSegment> CRAF_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION = dataLine -> {
        final String contig = dataLine.get(CRAFSegmentTableColumn.CONTIG);
        final int start = dataLine.getInt(CRAFSegmentTableColumn.START);
        final int end = dataLine.getInt(CRAFSegmentTableColumn.END);
        final int numPointsCopyRatio = dataLine.getInt(CRAFSegmentTableColumn.NUM_POINTS_COPY_RATIO);
        final int numPointsAlleleFraction = dataLine.getInt(CRAFSegmentTableColumn.NUM_POINTS_ALLELE_FRACTION);
        final double meanLog2CopyRatio = dataLine.getDouble(CRAFSegmentTableColumn.MEAN_LOG2_COPY_RATIO);
        final SimpleInterval interval = new SimpleInterval(contig, start, end);
        return new CRAFSegment(interval, numPointsCopyRatio, numPointsAlleleFraction, meanLog2CopyRatio);
    };

    private static final BiConsumer<CRAFSegment, DataLine> CRAF_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER = (alleleFractionSegment, dataLine) ->
            dataLine.append(alleleFractionSegment.getContig())
                    .append(alleleFractionSegment.getStart())
                    .append(alleleFractionSegment.getEnd())
                    .append(alleleFractionSegment.getNumPointsCopyRatio())
                    .append(alleleFractionSegment.getNumPointsAlleleFraction())
                    .append(alleleFractionSegment.getMeanLog2CopyRatio());

    public CRAFSegmentCollection(final File inputFile) {
        super(inputFile, CRAFSegmentTableColumn.COLUMNS, CRAF_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION, CRAF_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }

    public CRAFSegmentCollection(final String sampleName,
                                 final List<CRAFSegment> crafSegments) {
        super(sampleName, crafSegments, CRAFSegmentTableColumn.COLUMNS, CRAF_SEGMENT_DATA_LINE_TO_RECORD_FUNCTION, CRAF_SEGMENT_RECORD_AND_DATA_LINE_BI_CONSUMER);
    }

    public static CRAFSegmentCollection unionSegments(final CopyRatioSegmentCollection copyRatioSegments,
                                                      final CopyRatioCollection denoisedCopyRatios,
                                                      final AlleleFractionSegmentCollection alleleFractionSegments,
                                                      final AllelicCountCollection allelicCounts,
                                                      final int numCopyRatioIntervalsSmallSegmentThreshold) {
        Utils.validateArg(!(copyRatioSegments == null && alleleFractionSegments == null),
                "Must provide at least a copy-ratio segmentation or an allele-fraction segmentation.");
        ParamUtils.isPositiveOrZero(numCopyRatioIntervalsSmallSegmentThreshold,
                "Threshold number of copy-ratio intervals for small-segment merging must be non-negative.");

        final String sampleName;
        final List<CRAFSegment> crafSegments;

        if (alleleFractionSegments == null) {
            Utils.nonNull(denoisedCopyRatios);
            sampleName = copyRatioSegments.getSampleName();
            crafSegments = copyRatioSegments.getRecords().stream()
                    .map(s -> new CRAFSegment(s.getInterval(), s.getNumPoints(), 0, s.getMeanLog2CopyRatio()))
                    .collect(Collectors.toList());
        } else if (copyRatioSegments == null) {
            Utils.nonNull(allelicCounts);
            sampleName = alleleFractionSegments.getSampleName();
            crafSegments = alleleFractionSegments.getRecords().stream()
                    .map(s -> new CRAFSegment(s.getInterval(), 0, s.getNumPoints(), Double.NaN))
                    .collect(Collectors.toList());
        } else {
            Utils.validateArg(copyRatioSegments.getSampleName().equals(alleleFractionSegments.getSampleName()),
                    "Sample names from copy-ratio segmentation and allele-fraction segmentation must match.");
            sampleName = copyRatioSegments.getSampleName();

            //use old code for segment union
            logger.info(String.format("Combining %d copy-ratio segments and %d allele-fraction segments...",
                    copyRatioSegments.getRecords().size(), alleleFractionSegments.getRecords().size()));
            final Genome genome = new Genome(denoisedCopyRatios, allelicCounts);
            final List<SimpleInterval> unionedSegments = SegmentUtils.unionSegments(copyRatioSegments.getIntervals(), alleleFractionSegments.getIntervals(), genome);
            logger.info(String.format("After combining segments, %d segments remain...", unionedSegments.size()));
            logger.info(String.format("Merging segments with less than %d copy-ratio intervals...", numCopyRatioIntervalsSmallSegmentThreshold));
            final SegmentedGenome segmentedGenomeWithSmallSegments = new SegmentedGenome(unionedSegments, genome);
            final SegmentedGenome segmentedGenome = segmentedGenomeWithSmallSegments.mergeSmallSegments(numCopyRatioIntervalsSmallSegmentThreshold);
            logger.info(String.format("After merging small segments, %d segments remain...", segmentedGenome.getSegments().size()));

//            //union copy-ratio and allele-fraction segments
//            logger.info(String.format("Combining %d copy-ratio segments and %d allele-fraction segments...",
//                    copyRatioSegments.getRecords().size(), alleleFractionSegments.getRecords().size()));
//            crafSegments = new SegmentUnioner(copyRatioSegments, denoisedCopyRatios, alleleFractionSegments, allelicCounts)
//                    .constructUnionedCRAFSegments();
//            logger.info(String.format("After combining segments, %d segments remain...", crafSegments.size()));

            final OverlapDetector<CopyRatio> copyRatioOverlapDetector = OverlapDetector.create(
                    denoisedCopyRatios.getRecords().stream()
                            .map(cr -> new CopyRatio(cr.getMidpoint(), cr.getLog2CopyRatioValue())) //map copy-ratio intervals to their midpoints so that each will be uniquely contained in a single segment
                            .collect(Collectors.toList()));
            final OverlapDetector<AllelicCount> allelicCountOverlapDetector = OverlapDetector.create(allelicCounts.getRecords());
            crafSegments = segmentedGenome.getSegments().stream()
                    .map(s -> new CRAFSegment(
                            s,
                            copyRatioOverlapDetector.getOverlaps(s).stream()
                                    .map(CopyRatio::getLog2CopyRatioValue)
                                    .collect(Collectors.toList()),
                            new ArrayList<>(allelicCountOverlapDetector.getOverlaps(s))))
                    .collect(Collectors.toList());
        }
        return new CRAFSegmentCollection(sampleName, crafSegments);
    }

    public SegmentedGenome convertToSegmentedGenome(final CopyRatioCollection denoisedCopyRatios,
                                                    final AllelicCountCollection allelicCounts) {
        final Genome genome = new Genome(denoisedCopyRatios, allelicCounts);
        return new SegmentedGenome(getIntervals(), genome);
    }

    private static final class SegmentUnioner {
        final CopyRatioSegmentCollection copyRatioSegments;
        final CopyRatioCollection denoisedCopyRatios;
        final AlleleFractionSegmentCollection alleleFractionSegments;
        final AllelicCountCollection allelicCounts;
        final OverlapDetector<CopyRatio> copyRatioOverlapDetector;
        final OverlapDetector<AllelicCount> allelicCountOverlapDetector;

        private SegmentUnioner(final CopyRatioSegmentCollection copyRatioSegments,
                               final CopyRatioCollection denoisedCopyRatios,
                               final AlleleFractionSegmentCollection alleleFractionSegments,
                               final AllelicCountCollection allelicCounts) {
            this.copyRatioSegments = copyRatioSegments;
            this.denoisedCopyRatios = denoisedCopyRatios;
            this.alleleFractionSegments = alleleFractionSegments;
            this.allelicCounts = allelicCounts;
            copyRatioOverlapDetector = OverlapDetector.create(denoisedCopyRatios.getRecords());
            allelicCountOverlapDetector = OverlapDetector.create(allelicCounts.getRecords());
        }

        private List<CRAFSegment> constructUnionedCRAFSegments() {
            final List<SimpleInterval> unionedSegments = constructUnionedSegments(copyRatioSegments, copyRatioOverlapDetector, alleleFractionSegments, allelicCountOverlapDetector);
            final OverlapDetector<CopyRatio> copyRatioOverlapDetector = OverlapDetector.create(
                    denoisedCopyRatios.getRecords().stream()
                            .map(cr -> new CopyRatio(cr.getMidpoint(), cr.getLog2CopyRatioValue())) //map copy-ratio intervals to their midpoints so that each will be uniquely contained in a single segment
                            .collect(Collectors.toList()));
            return unionedSegments.stream()
                    .map(s -> new CRAFSegment(
                            s,
                            copyRatioOverlapDetector.getOverlaps(s).stream()
                                    .map(CopyRatio::getLog2CopyRatioValue)
                                    .collect(Collectors.toList()),
                            new ArrayList<>(allelicCountOverlapDetector.getOverlaps(s))))
                    .collect(Collectors.toList());
        }

        /**
         * Returns segments derived from the union of copy-ratio and allele-fraction segments.
         * First, all breakpoints from both sets of segments are combined to form new segments.
         * Spurious segments (i.e., segments containing only copy-ratio intervals that are created by
         * allele-fraction breakpoints and are not present in the original set of copy-ratio segments)
         * at the starts and ends of the original copy-ratio segments are then remerged to the right and left, respectively;
         * spurious segments introduced within the original copy-ratio segments are merged with adjacent segments
         * by removing the less favorable breakpoint according to {@link CopyRatioKernelSegmenter}.
         * Finally, the segments are trimmed by {@link #trimSegments}.
         */
        private static List<SimpleInterval> constructUnionedSegments(final CopyRatioSegmentCollection copyRatioSegments,
                                                                     final OverlapDetector<CopyRatio> copyRatioOverlapDetector,
                                                                     final AlleleFractionSegmentCollection alleleFractionSegments,
                                                                     final OverlapDetector<AllelicCount> allelicCountOverlapDetector) {
            final SortedMap<String, List<Breakpoint>> breakpointsByContig = collectBreakpointsByContig(copyRatioSegments, alleleFractionSegments);
            final List<SimpleInterval> untrimmedSegments = constructUntrimmedSegments(copyRatioOverlapDetector, allelicCountOverlapDetector, breakpointsByContig);

            //merge spurious untrimmed segments at copy-ratio-segment starts/ends created by allele-fraction breakpoints
            List<SimpleInterval> spuriousStartsAndEndsMergedSegments = mergeSpuriousStartsAndEnds(untrimmedSegments, copyRatioSegments, allelicCountOverlapDetector);
            //merge spurious untrimmed segments at allele-fraction-segment starts/ends created by copy-ratio breakpoints
            spuriousStartsAndEndsMergedSegments = mergeSpuriousStartsAndEnds(spuriousStartsAndEndsMergedSegments, alleleFractionSegments, copyRatioOverlapDetector);

            final List<SimpleInterval> spuriousMiddlesMergedSegments = mergeSpuriousMiddles(spuriousStartsAndEndsMergedSegments, copyRatioSegments, copyRatioOverlapDetector, allelicCountOverlapDetector);
            return spuriousMiddlesMergedSegments.stream().map(s -> trimSegments(s, copyRatioOverlapDetector, allelicCountOverlapDetector))
                    .collect(Collectors.toList());
        }

        private enum BreakpointType {
            START, END
        }

        private static final class Breakpoint extends Interval {
            final BreakpointType type;

            Breakpoint(final BreakpointType type, final String contig, final int site) {
                super(contig, site, site);
                this.type = type;
            }
            public BreakpointType getType() {
                return type;
            }

            public int getSite() {
                return getStart();
            }
        }

        /**
         * Returns a map of contig -> combined breakpoints from copy-ratio and allele-fraction segments on that contig.
         */
        private static SortedMap<String, List<Breakpoint>> collectBreakpointsByContig(final CopyRatioSegmentCollection copyRatioSegments,
                                                                                      final AlleleFractionSegmentCollection alleleFractionSegments) {
            return ListUtils.union(copyRatioSegments.getIntervals(), alleleFractionSegments.getIntervals()).stream()
                    .map(s -> Arrays.asList(
                            new Breakpoint(BreakpointType.START, s.getContig(), s.getStart()),
                            new Breakpoint(BreakpointType.END, s.getContig(), s.getEnd())))
                    .flatMap(Collection::stream)
                    .sorted()
                    .collect(Collectors.groupingBy(Breakpoint::getContig, TreeMap::new, Collectors.toList()));
        }

        /**
         * Given breakpointsByContig map, constructs a list of untrimmed segments (i.e., segments that are directly
         * adjacent to each other; segment boundaries at each breakpoint are determined by whether that breakpoint
         * was a start or an end for the corresponding original copy-ratio/allele-fraction segment),
         * then returns only those untrimmed segments that are non-empty.
         */
        private static List<SimpleInterval> constructUntrimmedSegments(final OverlapDetector<CopyRatio> copyRatioOverlapDetector,
                                                                       final OverlapDetector<AllelicCount> allelicCountOverlapDetector,
                                                                       final SortedMap<String, List<Breakpoint>> breakpointsByContig) {
            final List<SimpleInterval> segments = new ArrayList<>();
            for (final String contig : breakpointsByContig.keySet()) {
                final List<Breakpoint> breakpoints = breakpointsByContig.get(contig);
                int start = breakpoints.get(0).getSite();
                for (int i = 1; i < breakpoints.size(); i++) {
                    final Breakpoint breakpoint = breakpoints.get(i);
                    //adjust segment boundaries according to breakpoint type; this prevents, e.g., allelic-count sites
                    //that originally started allele-fraction segments from being stranded
                    final int end = breakpoint.getSite() - (breakpoint.type == BreakpointType.START ? 1 : 0);

                    if (end < start) {  //this could happen if there are adjacent breakpoints
                        continue;
                    }

                    final SimpleInterval segment = new SimpleInterval(contig, start, end);
                    if (copyRatioOverlapDetector.overlapsAny(segment) || allelicCountOverlapDetector.overlapsAny(segment)) {
                        segments.add(segment);
                    }
                    start = segment.getEnd() + 1;
                }
            }
            return segments;
        }

        /**
         * Given segments, returns a list with spurious segments---which contain only points of the type associated with
         * {@code originalPointSegments} and were created during segment union at the starts and ends of the
         * {@code originalPointSegments} by breakpoints arising from points of the other type---remerged.
         */
        private static <SEGMENT extends Locatable, POINT extends Locatable> List<SimpleInterval> mergeSpuriousStartsAndEnds(
                final List<SimpleInterval> segments,
                final TSVLocatableCollection<SEGMENT> originalPointSegments,
                final OverlapDetector<POINT> otherOverlapDetector) {

            //get starts and ends of original segments
            final Set<SimpleInterval> originalPointSegmentStarts =
                    originalPointSegments.getIntervals().stream()
                            .map(s -> new SimpleInterval(s.getContig(), s.getStart(), s.getStart()))
                            .collect(Collectors.toSet());
            final Set<SimpleInterval> originalPointSegmentEnds =
                    originalPointSegments.getIntervals().stream()
                            .map(s -> new SimpleInterval(s.getContig(), s.getEnd(), s.getEnd()))
                            .collect(Collectors.toSet());

            final List<SimpleInterval> mergedSegments = new ArrayList<>();
            final ListIterator<SimpleInterval> segmentsIter = segments.listIterator();
            while (segmentsIter.hasNext()) {
                final SimpleInterval segment = segmentsIter.next();
                //do not remerge segments containing points of the other type
                if (otherOverlapDetector.overlapsAny(segment)) {
                    mergedSegments.add(segment);
                    continue;
                }
                final SimpleInterval segmentStart =
                        new SimpleInterval(segment.getContig(), segment.getStart(), segment.getStart());
                final SimpleInterval segmentEnd =
                        new SimpleInterval(segment.getContig(), segment.getEnd(), segment.getEnd());
                if (originalPointSegmentStarts.contains(segmentStart) && !originalPointSegmentEnds.contains(segmentEnd)) {
                    //remerge segments introduced at starts to the right
                    final SimpleInterval nextSegment = segmentsIter.next();
                    mergedSegments.add(mergeSegments(segment, nextSegment));
                } else if (!originalPointSegmentStarts.contains(segmentStart) && originalPointSegmentEnds.contains(segmentEnd)) {
                    //remerge segments introduced at ends to the left
                    final int previousIndex = mergedSegments.size() - 1;
                    final SimpleInterval previousSegment = mergedSegments.get(previousIndex);
                    mergedSegments.set(previousIndex, mergeSegments(previousSegment, segment));
                } else {
                    //do not merge otherwise
                    mergedSegments.add(segment);
                }
            }
            return mergedSegments;
        }

        /**
         * Identifies spurious segments (i.e., segments containing only copy-ratio intervals that are created by
         * allele-fraction breakpoints and are not present in the original set of copy-ratio segments)
         * introduced into the middle of original copy-ratio segments by segment union
         * and merges them with adjacent segments by removing the less favorable breakpoint
         * according to {@link CopyRatioKernelSegmenter}, returning a new, modifiable list of segments.
         */
        private static List<SimpleInterval> mergeSpuriousMiddles(final List<SimpleInterval> segments,
                                                                 final CopyRatioSegmentCollection copyRatioSegments,
                                                                 final OverlapDetector<CopyRatio> copyRatioOverlapDetector,
                                                                 final OverlapDetector<AllelicCount> allelicCountOverlapDetector) {
//        final OverlapDetector<CopyRatio> copyRatioOverlapDetector = OverlapDetector.create(denoisedCopyRatios.getRecords());
//        final OverlapDetector<AllelicCount> allelicCountOverlapDetector = OverlapDetector.create(allelicCounts.getRecords());
//        final Set<SimpleInterval> copyRatioSegmentsSet = new HashSet<>(copyRatioSegments.getIntervals());
//        final List<SimpleInterval> mergedSegments = new ArrayList<>(segments);
//        int index = 0;
//        while (index < mergedSegments.size()) {
//            final SimpleInterval segment = mergedSegments.get(index);
//            final int numSNPs = genome.getSNPs().targetCount(segment);
//            //if current segment is a spurious middle, merge it with an adjacent segment
//            if (numSNPs == 0 && !copyRatioSegmentsSet.contains(segment)) {
//                final MergeDirection direction =
//                        SmallSegments.calculateMergeDirection(mergedSegments, genome, index);
//                if (direction == MergeDirection.LEFT) {
//                    //current = merge(left, current), remove left, stay on current during next iteration
//                    mergedSegments.set(index, mergeSegments(mergedSegments.get(index - 1), segment));
//                    mergedSegments.remove(index - 1);
//                    index -= 2;
//                } else if (direction == MergeDirection.RIGHT) {
//                    //current = merge(current, right), remove right, stay on current during next iteration
//                    mergedSegments.set(index, mergeSegments(segment, mergedSegments.get(index + 1)));
//                    mergedSegments.remove(index + 1);
//                    index--;
//                }
//            }
//            index++; //if no merge performed, go to next segment during next iteration
//        }
//        return mergedSegments;
            return segments;
        }

        /**
         * Given a segment and collections of copy ratios and allelic counts, returns a trimmed segment produced by
         * removing the empty portions at the start and the end of the original segment that do not overlap the
         * copy ratios and allelic counts that overlap with the original segment.
         * If this procedure would remove the entire segment, the original segment is returned instead.
         */
        private static SimpleInterval trimSegments(final SimpleInterval segment,
                                                   final OverlapDetector<CopyRatio> copyRatioOverlapDetector,
                                                   final OverlapDetector<AllelicCount> allelicCountOverlapDetector) {
//        final IndexRange targetRange = targets.indexRange(segment);
//        final IndexRange snpRange = snps.indexRange(segment);
//
//        final int numTargetsInInterval = targetRange.size();
//        final int numSNPsInInterval = snpRange.size();
//
//        int start = segment.getStart();
//        int end = segment.getEnd();
//
//        if (numTargetsInInterval == 0 && numSNPsInInterval > 0) {
//            //if there are no targets overlapping interval, use SNPs to determine trimmed interval
//            start = snps.target(snpRange.from).getStart();
//            end = snps.target(snpRange.to - 1).getEnd();
//        } else if (numTargetsInInterval > 0) {
//            //if interval start does not fall within first target, use start of first target as start of trimmed interval
//            start = Math.max(start, targets.target(targetRange.from).getStart());
//            //if interval end does not fall within last target, use end of last target as end of trimmed interval
//            end = Math.min(end, targets.target(targetRange.to - 1).getEnd());
//
//            if (numSNPsInInterval > 0) {
//                //if there are also SNPs within interval, check to see if they give a larger trimmed interval
//                start = Math.min(start, snps.target(snpRange.from).getStart());
//                end = Math.max(end, snps.target(snpRange.to - 1).getEnd());
//            }
//        }
//        if (start < segment.getStart() || end > segment.getEnd() || end < start) {
//            throw new GATKException.ShouldNeverReachHereException("Something went wrong in trimming interval.");
//        }
//        return new SimpleInterval(segment.getContig(), start, end);
            return segment;
        }

        private static SimpleInterval mergeSegments(final SimpleInterval segment1,
                                                    final SimpleInterval segment2) {
            Utils.validateArg(segment1.getContig().equals(segment2.getContig()),
                    String.format("Cannot join segments %s and %s on different chromosomes.", segment1.toString(), segment2.toString()));
            final int start = Math.min(segment1.getStart(), segment2.getStart());
            final int end = Math.max(segment1.getEnd(), segment2.getEnd());
            return new SimpleInterval(segment1.getContig(), start, end);
        }
    }
}