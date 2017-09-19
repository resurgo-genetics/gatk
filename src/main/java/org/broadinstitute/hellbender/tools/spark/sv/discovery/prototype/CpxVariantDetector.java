package org.broadinstitute.hellbender.tools.spark.sv.discovery.prototype;

import htsjdk.samtools.SAMSequenceDictionary;
import org.apache.logging.log4j.Logger;
import org.apache.spark.api.java.JavaRDD;
import org.apache.spark.broadcast.Broadcast;
import org.broadinstitute.hellbender.engine.datasources.ReferenceMultiSource;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignedContig;
import org.broadinstitute.hellbender.tools.spark.sv.discovery.AlignmentInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVInterval;
import org.broadinstitute.hellbender.tools.spark.sv.utils.SVLocation;

import java.util.ArrayList;
import java.util.List;
import java.util.ListIterator;
import java.util.stream.Collectors;

/**
 * This deals with the special case where a contig's multiple (> 2) alignments has head and tail mapped to the same chr.
 * For the case where the head and tail mapped to different chromosome, we could decide to emit all BND records, but
 * that could be dealt with later.
 */
final class CpxVariantDetector implements VariantDetectorFromLocalAssemblyContigAlignments {

    @Override
    public void inferSvAndWriteVCF(final JavaRDD<AlignedContig> localAssemblyContigs, final String vcfOutputFileName,
                                   final Broadcast<ReferenceMultiSource> broadcastReference, final SAMSequenceDictionary refSequenceDictionary,
                                   final Logger toolLogger){
        // preprocess AI list

        // after the preprocessing, the configuration might be different from the input configuration,
        // so again divert the contigs into different units

        // extract reference ordered jumping locations on reference

        // segment affect reference regions by jumping locations

        // make sense of event, i.e. provide interpretation, and extract corresponding alt haplotype

        // output VCF
    }


    private static List<AlignmentInterval> preprocessAlignments(final List<AlignmentInterval> originalConfiguration,
                                                                final int mapQThresholdInclusive,
                                                                final int uniqReadLenInclusive) {

        // two pass, each focusing on removing the alignments of a contig that offers low uniqueness

        // first pass is for removing alignments with low ref uniqueness, using low mapping quality as the criteria
        final List<AlignmentInterval> pickedAlignments =
                originalConfiguration.stream().filter(ai -> ai.mapQual >= mapQThresholdInclusive).collect(Collectors.toList());

        // second pass, the slower one, is to remove alignments offering low read uniqueness,
        // i.e. with only a very short part of the read being explained by this particular alignment;
        // search bi-directionally until cannot find overlap any more, subtract from it all overlaps.
        // This gives unique read region it explains. If this unique read region is "short" (e.g. less than half of
        // its aligner-assigned read consumption length, or shorter than 30 bp), drop it.
        final List<Integer> toRemove = new ArrayList<>(pickedAlignments.size());
        for (int i = 0; i < pickedAlignments.size(); ++i) { // implementation is less-than efficient because we are computing the overlap twice, but prototyping for now
            final AlignmentInterval cur = pickedAlignments.get(i);
            int maxOverlap = -1;
            for (int j = 0; j != i && j < pickedAlignments.size(); ++j) { // bi-directional
                final int overlap = AlignmentInterval.overlapOnContig(cur, pickedAlignments.get(j));
                if (overlap > 0)
                    maxOverlap = Math.max(maxOverlap, overlap);
                else // following ones, as guaranteed by the ordering of alignments in the contig, cannot overlap
                    break;
            }
            // TODO: 9/21/17 clip a copy of cur here and check if it is too short
            if (cur.endInAssembledContig - cur.startInAssembledContig + 1 < uniqReadLenInclusive)
                toRemove.add(i);
        }

        if ( toRemove.isEmpty() )
            return pickedAlignments;

        // removing in reverse order so that iterators are not invalidated if we were to remove from start
        final ListIterator<Integer> rit = toRemove.listIterator(toRemove.size());
        while (rit.hasPrevious()) {
            pickedAlignments.remove( rit.previous().intValue() );
        }

        return pickedAlignments;
    }

    /**
     * Each pair of neighboring reference locations are meant to be used closed, i.e. [a, b].
     */
    private static List<SVLocation> extractReferenceOrdereredJumpLocations(final List<AlignmentInterval> alignmentConfiguration) {

        // A configuration has a series jumps on the reference as indicated by the chimeric alignments.
        // A jump has a starting and landing ref location.
        // A jump can be: 1) gapped, meaning a part of read is uncovered by neighboring AI's;
        //                2) connected, meaning neighboring--but not overlapping on the read--AI's leave no base on the read uncovered;
        //                3) retracting, meaning neighboring AI's overlap on the read, pointing to homology between their ref span
        // Among them, retracting jumps are the most difficult to deal with, mainly due to how to have a consistent
        // homology-yielding scheme

        return null;
    }

    private static List<SVInterval> segmentReference(final List<SVLocation> jumpingLocations ) {
        return null;
    }


}
