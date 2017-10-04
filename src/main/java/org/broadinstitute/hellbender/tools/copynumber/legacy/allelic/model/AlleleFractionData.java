package org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.model;

import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCountCollection;
import org.broadinstitute.hellbender.tools.copynumber.legacy.formats.TSVLocatableCollection;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.DataCollection;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * {@link DataCollection} for the allele-fraction model containing the het alt and ref counts grouped by segment.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class AlleleFractionData implements DataCollection {
    private final AllelicCountCollection allelicCounts;
    private final List<SimpleInterval> segments;

    private final List<List<IndexedAllelicCount>> indexedAllelicCountsPerSegment = new ArrayList<>();

    AlleleFractionData(final AllelicCountCollection allelicCounts,
                       final List<SimpleInterval> segments) {
        this.allelicCounts = Utils.nonNull(allelicCounts);
        this.segments = Utils.nonEmpty(segments).stream().sorted(TSVLocatableCollection.LEXICOGRAPHICAL_ORDER_COMPARATOR).collect(Collectors.toList());

        final OverlapDetector<AllelicCount> allelicCountOverlapDetector = allelicCounts.getOverlapDetector();
        int index = 0;
        for (final SimpleInterval segment : segments) {
            final List<AllelicCount> allelicCountsInSegment = allelicCountOverlapDetector.getOverlaps(segment).stream()
                    .sorted(TSVLocatableCollection.LEXICOGRAPHICAL_ORDER_COMPARATOR)
                    .collect(Collectors.toList());
            final int segmentStartIndex = index;
            final List<IndexedAllelicCount> indexedAllelicCountsInSegment = IntStream.range(0, allelicCountsInSegment.size()).boxed()
                    .map(i -> new IndexedAllelicCount(allelicCountsInSegment.get(i), segmentStartIndex + i))
                    .collect(Collectors.toList());
            indexedAllelicCountsPerSegment.add(indexedAllelicCountsInSegment);
            index += allelicCountsInSegment.size();
        }
    }

    int getNumSegments() {
        return segments.size();
    }

    int getNumPoints() {
        return allelicCounts.size();
    }

    List<IndexedAllelicCount> getIndexedAllelicCountsInSegment(final int segment) {
        return Collections.unmodifiableList(indexedAllelicCountsPerSegment.get(segment));
    }

    static final class IndexedAllelicCount extends AllelicCount {
        private final int index;

        private IndexedAllelicCount(final AllelicCount allelicCount,
                                    final int index) {
            super(allelicCount.getInterval(), allelicCount.getRefReadCount(), allelicCount.getAltReadCount(), allelicCount.getRefNucleotide(), allelicCount.getAltNucleotide());
            this.index = index;
        }

        int getIndex() {
            return index;
        }
    }
}
