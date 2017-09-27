package org.broadinstitute.hellbender.tools.copynumber.legacy.allelic.segmentation;

import htsjdk.samtools.util.Locatable;
import org.apache.commons.math3.stat.inference.AlternativeHypothesis;
import org.apache.commons.math3.stat.inference.BinomialTest;
import org.broadinstitute.hellbender.tools.copynumber.allelic.alleliccount.AllelicCount;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionData;
import org.broadinstitute.hellbender.tools.exome.allelefraction.AlleleFractionInitializer;
import org.broadinstitute.hellbender.tools.exome.allelefraction.MinorAlleleFractionCache;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.param.ParamUtils;

import java.util.List;
import java.util.stream.Collectors;

public class AlleleFractionSegment implements Locatable {
    private static final double ERROR_RATE = 0.001;
    private static final double HOM_GENOTYPING_P_VALUE_THRESHOLD = 0.05;
    private static final double DEFAULT_ALLELIC_BIAS = 1.05;

    private final SimpleInterval interval;
    private final int numPoints;
    private final double meanMinorAlleleFraction;

    public AlleleFractionSegment(final SimpleInterval interval,
                                 final int numPoints,
                                 final double meanMinorAlleleFraction) {
        Utils.nonNull(interval);
        ParamUtils.isPositive(numPoints, "Number of points must be positive.");
        this.interval = interval;
        this.numPoints = numPoints;
        this.meanMinorAlleleFraction = meanMinorAlleleFraction;
    }

    public AlleleFractionSegment(final SimpleInterval interval,
                                 final List<AllelicCount> allelicCounts) {
        Utils.nonNull(interval);
        Utils.nonNull(allelicCounts);
        this.interval = interval;
        numPoints = allelicCounts.size();
        meanMinorAlleleFraction = allelicCounts.stream()
                .filter(ac -> new BinomialTest().binomialTest(
                        ac.getRefReadCount() + ac.getAltReadCount(),
                        Math.min(ac.getAltReadCount(), ac.getRefReadCount()),
                        ERROR_RATE,
                        AlternativeHypothesis.TWO_SIDED) <= HOM_GENOTYPING_P_VALUE_THRESHOLD)
                .mapToDouble(ac -> MinorAlleleFractionCache.get(ac.getAltReadCount(), ac.getRefReadCount(), DEFAULT_ALLELIC_BIAS))
                .average()
                .orElse(Double.NaN);
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

    public double getMeanMinorAlleleFraction() {
        return meanMinorAlleleFraction;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }
        final AlleleFractionSegment that = (AlleleFractionSegment) o;
        return numPoints == that.numPoints &&
                Double.compare(that.meanMinorAlleleFraction, meanMinorAlleleFraction) == 0 &&
                interval.equals(that.interval);
    }

    @Override
    public int hashCode() {
        int result;
        long temp;
        result = interval.hashCode();
        result = 31 * result + numPoints;
        temp = Double.doubleToLongBits(meanMinorAlleleFraction);
        result = 31 * result + (int) (temp ^ (temp >>> 32));
        return result;
    }
}
