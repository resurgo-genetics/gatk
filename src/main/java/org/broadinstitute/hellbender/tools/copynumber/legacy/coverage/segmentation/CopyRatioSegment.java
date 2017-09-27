package org.broadinstitute.hellbender.tools.copynumber.legacy.coverage.segmentation;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.List;

public class CopyRatioSegment implements Locatable {
    private final SimpleInterval interval;
    private final int numPoints;
    private final double meanLog2CopyRatio;

    public CopyRatioSegment(final SimpleInterval interval,
                            final int numPoints,
                            final double meanLog2CopyRatio) {
        Utils.nonNull(interval);
        ParamUtils.isPositive(numPoints, "Number of points must be positive.");
        this.interval = interval;
        this.numPoints = numPoints;
        this.meanLog2CopyRatio = meanLog2CopyRatio;
    }

    public CopyRatioSegment(final SimpleInterval interval,
                            final List<Double> denoisedCopyRatios) {
        this(interval, denoisedCopyRatios.size(), denoisedCopyRatios.stream().mapToDouble(Double::doubleValue).average().orElse(Double.NaN));
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

    public int getNumPoints() {
        return numPoints;
    }

    public double getMeanLog2CopyRatio() {
        return meanLog2CopyRatio;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final CopyRatioSegment that = (CopyRatioSegment) o;
        if (numPoints != that.numPoints) {
            return false;
        }
        if (Double.compare(that.meanLog2CopyRatio, meanLog2CopyRatio) != 0) {
            return false;
        }
        return interval.equals(that.interval);
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = interval.hashCode();
        result = 31 * result + numPoints;
        temp = Double.doubleToLongBits(meanLog2CopyRatio);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }
}
