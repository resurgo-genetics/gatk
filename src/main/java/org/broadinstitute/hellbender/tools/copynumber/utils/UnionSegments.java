package org.broadinstitute.hellbender.tools.copynumber.utils;

import com.google.common.collect.Lists;
import com.google.common.collect.Sets;
import htsjdk.samtools.SAMSequenceDictionary;
import htsjdk.samtools.util.Locatable;
import org.broadinstitute.barclay.argparser.Argument;
import org.broadinstitute.barclay.argparser.CommandLineProgramProperties;
import org.broadinstitute.hellbender.cmdline.ExomeStandardArgumentDefinitions;
import org.broadinstitute.hellbender.cmdline.programgroups.CopyNumberProgramGroup;
import org.broadinstitute.hellbender.engine.GATKTool;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.SimpleAnnotatedGenomicRegion;
import org.broadinstitute.hellbender.tools.copynumber.utils.annotatedregion.VersatileAnnotatedRegionParser;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;

import java.io.File;
import java.io.IOException;
import java.util.*;
import java.util.function.Function;
import java.util.stream.Collectors;

@CommandLineProgramProperties(
        oneLineSummary = "Do a breakpoint union of two segment files and annotate with chosen columns from each file.",
        summary = "Breakpoint union of two segment files while preserving annotations.\n" +
                "This tool will load all segments into RAM.",
        programGroup = CopyNumberProgramGroup.class)
public class UnionSegments extends GATKTool {

    @Argument(
            doc = "Input segment files -- must be specified twice, but order does not matter.",
            fullName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_LONG_NAME,
            shortName = ExomeStandardArgumentDefinitions.SEGMENT_FILE_SHORT_NAME,
            maxElements = 2, minElements = 2
    )
    protected List<File> segmentsFile;

    @Argument(
            doc="List of columns in either segment file that should be reported in the output file.  If the column header exists in both, it will have an append.",
            fullName = "columnsOfInterest", shortName = "cols", minElements = 1
    )
    protected Set<String> columnsOfInterest;

    @Override
    public boolean requiresReference() {
        return true;
    }

    @Override
    public void traverse() {
        try {
            final VersatileAnnotatedRegionParser parser = new VersatileAnnotatedRegionParser();
            final List<SimpleAnnotatedGenomicRegion> segments1 = parser.readAnnotatedRegions(segmentsFile.get(0), columnsOfInterest);
            final List<SimpleAnnotatedGenomicRegion> segments2 = parser.readAnnotatedRegions(segmentsFile.get(1), columnsOfInterest);

            // Create a map of input to output headers.  I.e. the annotation in the segment1 to the output that should be written in the final file.
            //  This assumes that the keys in each entry of the list is the same.
            // TODO: If we want to support more than two segment files, this is the only bit that requires thinking.
            final Set<String> intersectingAnnotations = Sets.intersection(segments1.get(0).getAnnotations().keySet(), segments2.get(0).getAnnotations().keySet());

            // Create the obvious mappings that are identity then tack on new annotations for conflicts.
            final Map<String, String> input1ToOutputHeaderMap = Utils.stream(segments1.get(0).getAnnotations().keySet().iterator()).filter(a -> !intersectingAnnotations.contains(a))
                    .collect(Collectors.toMap(Function.identity(), Function.identity()));
            Utils.stream(intersectingAnnotations.iterator()).forEach(a -> input1ToOutputHeaderMap.put(a, a+"_1"));

            final Map<String, String> input2ToOutputHeaderMap = Utils.stream(segments2.get(0).getAnnotations().keySet().iterator()).filter(a -> !intersectingAnnotations.contains(a))
                    .collect(Collectors.toMap(Function.identity(), Function.identity()));
            Utils.stream(intersectingAnnotations.iterator()).forEach(a -> input2ToOutputHeaderMap.put(a, a+"_2"));

            final List<SimpleAnnotatedGenomicRegion> finalList = annotateUnionedIntervals(segments1, segments2,
                    getBestAvailableSequenceDictionary(), Lists.newArrayList(input1ToOutputHeaderMap, input2ToOutputHeaderMap));

            // TODO: Write output

        } catch (final IOException ioe) {
            throw new UserException.BadInput("Could not parse input file", ioe);
        }
    }

    private List<SimpleAnnotatedGenomicRegion> annotateUnionedIntervals(final List<SimpleAnnotatedGenomicRegion> segments1, final List<SimpleAnnotatedGenomicRegion> segments2,
                                                                        final SAMSequenceDictionary dictionary, final List<Map<String, String>> inputToOutputHeaderMaps) {

        // We assume that the union'ed intervals are sorted.
        final List<Locatable> unionIntervals = IntervalUtils.unionIntervals(segments1.stream().map(s -> s.getInterval()).collect(Collectors.toList()),
                segments2.stream().map(s -> s.getInterval()).collect(Collectors.toList()));

        final List<List<SimpleAnnotatedGenomicRegion>> segmentList = Lists.newArrayList(segments1, segments2);
        final List<Map<Locatable, List<SimpleAnnotatedGenomicRegion>>> unionIntervalsToSegmentsMaps = segmentList.stream()
                .map(segs -> IntervalUtils.createOverlapMap(unionIntervals, segs, dictionary)).collect(Collectors.toList());

        final List<SimpleAnnotatedGenomicRegion> result = new ArrayList<>();

        for (final Locatable interval: unionIntervals) {
            final Map<String, String> intervalAnnotationMap = new HashMap<>();

            for (int i = 0; i< unionIntervalsToSegmentsMaps.size(); i ++) {
                final Map<Locatable, List<SimpleAnnotatedGenomicRegion>> unionIntervalsToSegmentsMap = unionIntervalsToSegmentsMaps.get(i);
                final SimpleAnnotatedGenomicRegion seg = unionIntervalsToSegmentsMap.get(interval).get(0);
                inputToOutputHeaderMaps.get(i).entrySet().stream()
                        .forEach(e -> intervalAnnotationMap.put(e.getValue(), seg.getAnnotations().get(e.getKey())));
            }

            result.add(new SimpleAnnotatedGenomicRegion(new SimpleInterval(interval), intervalAnnotationMap));
        }

        return result;
    }



}
