package org.broadinstitute.hellbender.tools.copynumber.legacy.formats;

import htsjdk.samtools.util.Locatable;
import htsjdk.samtools.util.OverlapDetector;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.tsv.DataLine;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.Comparator;
import java.util.List;
import java.util.function.BiConsumer;
import java.util.function.Function;
import java.util.stream.Collectors;

/**
 * Represents a coordinate-sorted collection of records that extend {@link Locatable}
 * (although contigs are assumed to be non-null when writing to file) associated with a sample name, a set of
 * mandatory column headers given by a {@link TableColumnCollection}, and lambdas for
 * reading and writing records.  Records are sorted using {@link IntervalUtils#LEXICOGRAPHICAL_ORDER_COMPARATOR}.
 * See TSVLocatableCollectionUnitTest for a simple example of a subclass.
 *
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public abstract class LocatableCollectionTSV<T extends Locatable> extends NamedSampleTSV<T> {
    public static final Comparator<Locatable> LEXICOGRAPHICAL_ORDER_COMPARATOR = IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR;

    /**
     * Records are sorted using {@code LEXICOGRAPHICAL_ORDER_COMPARATOR}.
     */
    public LocatableCollectionTSV(final String sampleName,
                                  final List<T> records,
                                  final TableColumnCollection mandatoryColumns,
                                  final Function<DataLine, T> dataLineToRecordFunction,
                                  final BiConsumer<T, DataLine> recordAndDataLineBiConsumer) {
        super(
                Utils.nonNull(sampleName),
                Utils.nonNull(records).stream().sorted(LEXICOGRAPHICAL_ORDER_COMPARATOR).collect(Collectors.toList()),
                Utils.nonNull(mandatoryColumns),
                Utils.nonNull(dataLineToRecordFunction),
                Utils.nonNull(recordAndDataLineBiConsumer));
    }

    /**
     * @throws UserException.BadInput if records are not sorted using {@code LEXICOGRAPHICAL_ORDER_COMPARATOR}
     */
    public LocatableCollectionTSV(final File inputFile,
                                  final TableColumnCollection mandatoryColumns,
                                  final Function<DataLine, T> dataLineToRecordFunction,
                                  final BiConsumer<T, DataLine> recordAndDataLineBiConsumer) {
        super(
                Utils.nonNull(inputFile),
                Utils.nonNull(mandatoryColumns),
                Utils.nonNull(dataLineToRecordFunction),
                Utils.nonNull(recordAndDataLineBiConsumer));
        if (!getRecords().stream().sorted(LEXICOGRAPHICAL_ORDER_COMPARATOR).collect(Collectors.toList()).equals(getRecords())) {
            throw new UserException.BadInput(String.format("Records in input file %s were not sorted in lexicographical order.", inputFile));
        }
    }

    /**
     * @return  a new modifiable list of {@link SimpleInterval}s corresponding to the {@link Locatable}s
     *          for each record contained in the collection
     */
    public List<SimpleInterval> getIntervals() {
        return getRecords().stream()
                .map(r -> new SimpleInterval(r.getContig(), r.getStart(), r.getEnd()))
                .collect(Collectors.toList());
    }

    public OverlapDetector<T> getOverlapDetector() {
        return OverlapDetector.create(getRecords());
    }
}