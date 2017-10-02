package org.broadinstitute.hellbender.tools.exome.allelefraction;

import org.broadinstitute.hellbender.tools.exome.SegmentedGenome;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.mcmc.DataCollection;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

/**
 * {@link DataCollection} for the allele-fraction model containing the set of het alt and ref counts
 * and the grouping of hets into segments.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionData implements DataCollection {
    private final int numSegments;
    private final int numPoints;
    private final AllelicPanelOfNormals allelicPoN;
    private final List<List<AllelicCount>> allelicCountsPerSegment;

    public AlleleFractionData(final SegmentedGenome segmentedGenome) {
        this(segmentedGenome, AllelicPanelOfNormals.EMPTY_PON);
    }

    public AlleleFractionData(final SegmentedGenome segmentedGenome, final AllelicPanelOfNormals allelicPoN) {
        numSegments = segmentedGenome.getSegments().size();
        this.allelicPoN = allelicPoN;
        allelicCountsPerSegment = new ArrayList<>(segmentedGenome.getSegments().size());
        final TargetCollection<AllelicCount> alleleCounts = segmentedGenome.getGenome().getSNPs();
        int numPoints = 0;
        for (final SimpleInterval segment : segmentedGenome.getSegments()) {
            final List<AllelicCount> allelicCountsInSegment = alleleCounts.targets(segment);
            allelicCountsPerSegment.add(allelicCountsInSegment);
            numPoints += allelicCountsInSegment.size();
        }
        this.numPoints = numPoints;
    }

    public AllelicPanelOfNormals getPoN() { return allelicPoN; }

    public List<AllelicCount> getCountsInSegment(final int segment) {
        return Collections.unmodifiableList(allelicCountsPerSegment.get(segment));
    }

    public int getNumHetsInSegment(final int segment) {
        return allelicCountsPerSegment.get(segment).size();
    }

    public int getNumSegments() { return numSegments; }

    public int getNumPoints() { return numPoints; }
}
