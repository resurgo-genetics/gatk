package org.broadinstitute.hellbender.utils.smithwaterman;

import htsjdk.samtools.Cigar;
import htsjdk.samtools.CigarElement;
import htsjdk.samtools.CigarOperator;
import htsjdk.samtools.TextCigarCodec;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWAlignerNativeBinding;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWNativeAlignerResult;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWOverhangStrategy;
import org.broadinstitute.gatk.nativebindings.smithwaterman.SWParameters;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.read.AlignmentUtils;

import java.util.ArrayList;
import java.util.List;

/**
 * A wrapper that converts instances of {@link SWAlignerNativeBinding} into a {@link SmithWatermanAligner}
 */
public final class SWNativeAlignerWrapper implements SmithWatermanAligner {
    private final SWAlignerNativeBinding aligner;

    public SWNativeAlignerWrapper(final SWAlignerNativeBinding aligner) {
        this.aligner = aligner;
    }

    /**
     * The state of a trace step through the matrix
     */
    protected enum State {
        MATCH,
        INSERTION,
        DELETION,
        CLIP
    }


    @Override
    public SmithWatermanAlignment align(final byte[] reference, final byte[] alternate, final SWParameters parameters, final SWOverhangStrategy overhangStrategy){

        Utils.nonNull(parameters);
        Utils.nonNull(overhangStrategy);

        // avoid running full Smith-Waterman if there is an exact match of alternate in reference
        int matchIndex = -1;
        if (overhangStrategy == SWOverhangStrategy.SOFTCLIP || overhangStrategy == SWOverhangStrategy.IGNORE) {
            // Use a substring search to find an exact match of the alternate in the reference
            // NOTE: This approach only works for SOFTCLIP and IGNORE overhang strategies
            matchIndex = Utils.lastIndexOf(reference, alternate);
        }

        if (matchIndex != -1) {
            // generate the alignment result when the substring search was successful
            final List<CigarElement> lce = new ArrayList<>(alternate.length);
            CigarOperator op = null;
            switch (State.MATCH) {
                case MATCH: op = CigarOperator.M; break;
                case INSERTION: op = CigarOperator.I; break;
                case DELETION: op = CigarOperator.D; break;
                case CLIP: op = CigarOperator.S; break;
            }
            lce.add(new CigarElement(alternate.length, op));
            return new SWNativeResultWrapper(AlignmentUtils.consolidateCigar(new Cigar(lce)), matchIndex);
        }
        else {
            // run full Smith-Waterman
            final int n = reference.length+1;
            final int m = alternate.length+1;
            final int[][] sw = new int[n][m];
            final int[][] btrack=new int[n][m];
            final SWNativeAlignerResult alignment = aligner.align(reference, alternate,parameters,overhangStrategy);

            return new SWNativeResultWrapper(alignment);
        }

    }

    @Override
    public void close() {
        aligner.close();
    }

    private static final class SWNativeResultWrapper implements SmithWatermanAlignment {
        private final Cigar cigar;
        private final int alignmentOffset;

        public SWNativeResultWrapper(final SWNativeAlignerResult nativeResult) {
            this.cigar = TextCigarCodec.decode(nativeResult.cigar);
            this.alignmentOffset = nativeResult.alignment_offset;
        }

        public SWNativeResultWrapper(final Cigar cigar, final int alignmentOffset) {
            this.cigar = cigar;
            this.alignmentOffset =  alignmentOffset;
        }

        @Override
        public Cigar getCigar() {
            return cigar;
        }

        @Override
        public int getAlignmentOffset() {
            return alignmentOffset;
        }
    }


}
