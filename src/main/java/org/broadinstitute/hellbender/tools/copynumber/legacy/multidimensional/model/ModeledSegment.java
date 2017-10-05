package org.broadinstitute.hellbender.tools.copynumber.legacy.multidimensional.model;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.Decile;
import org.broadinstitute.hellbender.utils.mcmc.PosteriorSummary;

public class ModeledSegment implements Locatable {
    private final SimpleInterval interval;
    private final int numPointsCopyRatio;
    private final int numPointsAlleleFraction;

    private final SimplePosteriorSummary log2CopyRatioSimplePosteriorSummary;
    private final SimplePosteriorSummary minorAlleleFractionSimplePosteriorSummary;

    public ModeledSegment(final SimpleInterval interval,
                          final int numPointsCopyRatio,
                          final int numPointsAlleleFraction,
                          final SimplePosteriorSummary log2CopyRatioSimplePosteriorSummary,
                          final SimplePosteriorSummary minorAlleleFractionSimplePosteriorSummary) {
        Utils.validateArg(numPointsCopyRatio > 0 || numPointsAlleleFraction > 0,
                "Number of copy-ratio points or number of allele-fraction points must be positive.");
        this.interval = Utils.nonNull(interval);
        this.numPointsCopyRatio = numPointsCopyRatio;
        this.numPointsAlleleFraction = numPointsAlleleFraction;
        this.log2CopyRatioSimplePosteriorSummary = Utils.nonNull(log2CopyRatioSimplePosteriorSummary);
        this.minorAlleleFractionSimplePosteriorSummary = Utils.nonNull(minorAlleleFractionSimplePosteriorSummary);
    }

    public ModeledSegment(final SimpleInterval interval,
                          final int numPointsCopyRatio,
                          final int numPointsAlleleFraction,
                          final PosteriorSummary log2CopyRatioPosteriorSummary,
                          final PosteriorSummary minorAlleleFractionPosteriorSummary) {
        this.interval = Utils.nonNull(interval);
        Utils.validateArg(numPointsCopyRatio > 0 || numPointsAlleleFraction > 0,
                "Number of copy-ratio points or number of allele-fraction points must be positive.");
        Utils.nonNull(log2CopyRatioPosteriorSummary);
        Utils.nonNull(minorAlleleFractionPosteriorSummary);
        this.numPointsCopyRatio = numPointsCopyRatio;
        this.numPointsAlleleFraction = numPointsAlleleFraction;
        log2CopyRatioSimplePosteriorSummary = new SimplePosteriorSummary(
                log2CopyRatioPosteriorSummary.getCenter(),
                log2CopyRatioPosteriorSummary.getDeciles().get(Decile.DECILE_10),
                log2CopyRatioPosteriorSummary.getDeciles().get(Decile.DECILE_50),
                log2CopyRatioPosteriorSummary.getDeciles().get(Decile.DECILE_90));
        minorAlleleFractionSimplePosteriorSummary = new SimplePosteriorSummary(
                minorAlleleFractionPosteriorSummary.getCenter(),
                minorAlleleFractionPosteriorSummary.getDeciles().get(Decile.DECILE_10),
                minorAlleleFractionPosteriorSummary.getDeciles().get(Decile.DECILE_50),
                minorAlleleFractionPosteriorSummary.getDeciles().get(Decile.DECILE_90));
    }

    @Override
    public String getContig() {
        return interval.getContig();
    }

    @Override
    public int getStart() {
        return interval.getStart();
    }

    @Override
    public int getEnd() {
        return interval.getEnd();
    }

    public SimpleInterval getInterval() {
        return interval;
    }

    public int getNumPointsCopyRatio() {
        return numPointsCopyRatio;
    }

    public int getNumPointsAlleleFraction() {
        return numPointsAlleleFraction;
    }

    public SimplePosteriorSummary getLog2CopyRatioSimplePosteriorSummary() {
        return log2CopyRatioSimplePosteriorSummary;
    }

    public SimplePosteriorSummary getMinorAlleleFractionSimplePosteriorSummary() {
        return minorAlleleFractionSimplePosteriorSummary;
    }

    public static final class SimplePosteriorSummary {
        private final double mode;
        private final double decile10;
        private final double decile50;
        private final double decile90;

        public SimplePosteriorSummary(final double mode,
                                      final double decile10,
                                      final double decile50,
                                      final double decile90) {
            this.mode = mode;
            this.decile10 = decile10;
            this.decile50 = decile50;
            this.decile90 = decile90;
        }

        public double getMode() {
            return mode;
        }

        public double getDecile10() {
            return decile10;
        }

        public double getDecile50() {
            return decile50;
        }

        public double getDecile90() {
            return decile90;
        }

        @Override
        public String toString() {
            return "SimplePosteriorSummary{" +
                    "mode=" + mode +
                    ", decile10=" + decile10 +
                    ", decile50=" + decile50 +
                    ", decile90=" + decile90 +
                    '}';
        }
    }
}
