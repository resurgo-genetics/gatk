package org.broadinstitute.hellbender.tools.copynumber.legacy.multidimensional.model;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.tools.exome.ACNVModeledSegment;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.Decile;

public class ModeledSegment implements Locatable {
    private final SimpleInterval interval;
    private final int numPointsCopyRatio;
    private final int numPointsAlleleFraction;
    private final double log2CopyRatioPosteriorMode;
    private final double log2CopyRatioPosterior10;
    private final double log2CopyRatioPosterior50;
    private final double log2CopyRatioPosterior90;

    private final double minorAlleleFractionPosteriorMode;
    private final double minorAlleleFractionPosterior10;
    private final double minorAlleleFractionPosterior50;
    private final double minorAlleleFractionPosterior90;

    public ModeledSegment(final SimpleInterval interval,
                          final int numPointsCopyRatio,
                          final int numPointsAlleleFraction,
                          final double log2CopyRatioPosteriorMode,
                          final double log2CopyRatioPosterior10,
                          final double log2CopyRatioPosterior50,
                          final double log2CopyRatioPosterior90,
                          final double minorAlleleFractionPosteriorMode,
                          final double minorAlleleFractionPosterior10,
                          final double minorAlleleFractionPosterior50,
                          final double minorAlleleFractionPosterior90) {
        Utils.nonNull(interval);
        Utils.validateArg(numPointsCopyRatio > 0 || numPointsAlleleFraction > 0,
                "Number of copy-ratio points or number of allele-fraction points must be positive.");
        this.interval = interval;
        this.numPointsCopyRatio = numPointsCopyRatio;
        this.numPointsAlleleFraction = numPointsAlleleFraction;
        this.log2CopyRatioPosteriorMode = log2CopyRatioPosteriorMode;
        this.log2CopyRatioPosterior10 = log2CopyRatioPosterior10;
        this.log2CopyRatioPosterior50 = log2CopyRatioPosterior50;
        this.log2CopyRatioPosterior90 = log2CopyRatioPosterior90;
        this.minorAlleleFractionPosteriorMode = minorAlleleFractionPosteriorMode;
        this.minorAlleleFractionPosterior10 = minorAlleleFractionPosterior10;
        this.minorAlleleFractionPosterior50 = minorAlleleFractionPosterior50;
        this.minorAlleleFractionPosterior90 = minorAlleleFractionPosterior90;
    }

    public ModeledSegment(final int numPointsCopyRatio,
                          final int numPointsAlleleFraction,
                          final ACNVModeledSegment acnvModeledSegment) {
        Utils.validateArg(numPointsCopyRatio > 0 || numPointsAlleleFraction > 0,
                "Number of copy-ratio points or number of allele-fraction points must be positive.");
        Utils.nonNull(acnvModeledSegment);
        this.interval = acnvModeledSegment.getInterval();
        this.numPointsCopyRatio = numPointsCopyRatio;
        this.numPointsAlleleFraction = numPointsAlleleFraction;
        this.log2CopyRatioPosteriorMode = acnvModeledSegment.getSegmentMeanPosteriorSummary().getCenter();
        this.log2CopyRatioPosterior10 = acnvModeledSegment.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_10);
        this.log2CopyRatioPosterior50 = acnvModeledSegment.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_50);
        this.log2CopyRatioPosterior90 = acnvModeledSegment.getSegmentMeanPosteriorSummary().getDeciles().get(Decile.DECILE_90);
        this.minorAlleleFractionPosteriorMode = acnvModeledSegment.getMinorAlleleFractionPosteriorSummary().getCenter();
        this.minorAlleleFractionPosterior10 = acnvModeledSegment.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_10);
        this.minorAlleleFractionPosterior50 = acnvModeledSegment.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_50);
        this.minorAlleleFractionPosterior90 = acnvModeledSegment.getMinorAlleleFractionPosteriorSummary().getDeciles().get(Decile.DECILE_90);
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

    public double getLog2CopyRatioPosteriorMode() {
        return log2CopyRatioPosteriorMode;
    }

    public double getLog2CopyRatioPosterior10() {
        return log2CopyRatioPosterior10;
    }

    public double getLog2CopyRatioPosterior50() {
        return log2CopyRatioPosterior50;
    }

    public double getLog2CopyRatioPosterior90() {
        return log2CopyRatioPosterior90;
    }

    public double getMinorAlleleFractionPosteriorMode() {
        return minorAlleleFractionPosteriorMode;
    }

    public double getMinorAlleleFractionPosterior10() {
        return minorAlleleFractionPosterior10;
    }

    public double getMinorAlleleFractionPosterior50() {
        return minorAlleleFractionPosterior50;
    }

    public double getMinorAlleleFractionPosterior90() {
        return minorAlleleFractionPosterior90;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final ModeledSegment that = (ModeledSegment) o;
        return numPointsCopyRatio == that.numPointsCopyRatio &&
                numPointsAlleleFraction == that.numPointsAlleleFraction &&
                Double.compare(that.log2CopyRatioPosteriorMode, log2CopyRatioPosteriorMode) == 0 &&
                Double.compare(that.log2CopyRatioPosterior10, log2CopyRatioPosterior10) == 0 &&
                Double.compare(that.log2CopyRatioPosterior50, log2CopyRatioPosterior50) == 0 &&
                Double.compare(that.log2CopyRatioPosterior90, log2CopyRatioPosterior90) == 0 &&
                Double.compare(that.minorAlleleFractionPosteriorMode, minorAlleleFractionPosteriorMode) == 0 &&
                Double.compare(that.minorAlleleFractionPosterior10, minorAlleleFractionPosterior10) == 0 &&
                Double.compare(that.minorAlleleFractionPosterior50, minorAlleleFractionPosterior50) == 0 &&
                Double.compare(that.minorAlleleFractionPosterior90, minorAlleleFractionPosterior90) == 0 &&
                interval.equals(that.interval);
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = interval.hashCode();
        result = 31 * result + numPointsCopyRatio;
        result = 31 * result + numPointsAlleleFraction;
        temp = Double.doubleToLongBits(log2CopyRatioPosteriorMode);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(log2CopyRatioPosterior10);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(log2CopyRatioPosterior50);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(log2CopyRatioPosterior90);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(minorAlleleFractionPosteriorMode);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(minorAlleleFractionPosterior10);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(minorAlleleFractionPosterior50);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        temp = Double.doubleToLongBits(minorAlleleFractionPosterior90);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }
}
