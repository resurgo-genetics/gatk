package org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.model;

import org.apache.commons.math3.random.RandomGenerator;
import org.apache.logging.log4j.LogManager;
import org.apache.logging.log4j.Logger;
import org.broadinstitute.hellbender.tools.exome.Genome;
import org.broadinstitute.hellbender.tools.exome.SegmentedGenome;
import org.broadinstitute.hellbender.tools.exome.TargetCollection;
import org.broadinstitute.hellbender.tools.exome.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.pon.allelic.AllelicPanelOfNormals;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.DataCollection;

import java.util.*;
import java.util.stream.Collectors;
import java.util.stream.IntStream;

/**
 * {@link DataCollection} for the allele-fraction model containing the set of het alt and ref counts
 * and the grouping of hets into segments.
 *
 * @author David Benjamin &lt;davidben@broadinstitute.org&gt;
 */
public final class AlleleFractionData implements DataCollection {
    private static final Logger logger = LogManager.getLogger(AlleleFractionData.class);

    private final SegmentedGenome segmentedGenome;
    private final int numSegments;
    private final int numPoints;
    private final AllelicPanelOfNormals allelicPoN;
    private final List<List<AllelicCount>> allelicCountsPerSegment;
    private final Map<Integer, IndexPair> indexToIndexPairMap;

    private static final class IndexPair {
        private final int segmentIndex;
        private final int countIndexWithinSegment;

        private IndexPair(final int segmentIndex, final int countIndexWithinSegment) {
            this.segmentIndex = segmentIndex;
            this.countIndexWithinSegment = countIndexWithinSegment;
        }
    }

    public AlleleFractionData(final SegmentedGenome segmentedGenome) {
        this(segmentedGenome, AllelicPanelOfNormals.EMPTY_PON);
    }

    public AlleleFractionData(final SegmentedGenome segmentedGenome, final AllelicPanelOfNormals allelicPoN) {
        this.segmentedGenome = segmentedGenome;
        this.allelicPoN = allelicPoN;
        numSegments = segmentedGenome.getSegments().size();
        indexToIndexPairMap = new HashMap<>();
        allelicCountsPerSegment = new ArrayList<>(segmentedGenome.getSegments().size());
        final TargetCollection<AllelicCount> alleleCounts = segmentedGenome.getGenome().getSNPs();
        int numPoints = 0;
        for (int segmentIndex = 0; segmentIndex < numSegments; segmentIndex++) {
            final int currentNumPoints = numPoints;
            final SimpleInterval segment = segmentedGenome.getSegments().get(segmentIndex);
            final List<AllelicCount> allelicCountsInSegment = alleleCounts.targets(segment);
            allelicCountsPerSegment.add(allelicCountsInSegment);
            final int numPointsInSegment = allelicCountsInSegment.size();
            for (int countIndexWithinSegment = 0; countIndexWithinSegment < numPointsInSegment; countIndexWithinSegment++) {
                indexToIndexPairMap.put(currentNumPoints + countIndexWithinSegment, new IndexPair(segmentIndex, countIndexWithinSegment));
            }
            numPoints += numPointsInSegment;
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

    AlleleFractionData subsample(final RandomGenerator rng, final int numPointsSubsampling) {
        Utils.nonNull(rng);
        Utils.validateArg(numPointsSubsampling <= numPoints, "Number of points to subsample must be less than number of points.");
        logger.info(String.format("Subsampling %d / %d points...", numPointsSubsampling, numPoints));
        final List<Integer> shuffledIndices = IntStream.range(0, numPoints).boxed().collect(Collectors.toList());
        Collections.shuffle(shuffledIndices, new Random(rng.nextInt()));
        final List<AllelicCount> subsampledCounts = shuffledIndices.subList(0, numPointsSubsampling).stream()
                .map(indexToIndexPairMap::get)
                .map(ip -> allelicCountsPerSegment.get(ip.segmentIndex).get(ip.countIndexWithinSegment))
                .collect(Collectors.toList());
        final Genome subsampledGenome = new Genome(segmentedGenome.getGenome().getTargets().targets(), subsampledCounts, segmentedGenome.getGenome().getSampleName());
        return new AlleleFractionData(new SegmentedGenome(segmentedGenome.getSegments(), subsampledGenome), allelicPoN);
    }
}
