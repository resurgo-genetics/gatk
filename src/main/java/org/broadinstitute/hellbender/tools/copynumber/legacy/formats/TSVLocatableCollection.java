package org.broadinstitute.hellbender.tools.copynumber.legacy.formats;

import htsjdk.samtools.util.Locatable;
import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.IntervalUtils;
import org.broadinstitute.hellbender.utils.SimpleInterval;
import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.tsv.*;

import java.io.File;
import java.io.IOException;
import java.io.UncheckedIOException;
import java.util.Collections;
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
public abstract class TSVLocatableCollection<T extends Locatable> {
    private final String sampleName;
    private final List<T> records;
    private final TableColumnCollection mandatoryColumns;
    private final Function<DataLine, T> dataLineToRecordFunction;
    private final BiConsumer<T, DataLine> recordAndDataLineBiConsumer;

    /**
     * Constructor given the sample name, the list of records, the mandatory column headers,
     * and the lambdas for reading and writing records.
     * Records are sorted using {@link IntervalUtils#LEXICOGRAPHICAL_ORDER_COMPARATOR}.
     */
    public TSVLocatableCollection(final String sampleName,
                                  final List<T> records,
                                  final TableColumnCollection mandatoryColumns,
                                  final Function<DataLine, T> dataLineToRecordFunction,
                                  final BiConsumer<T, DataLine> recordAndDataLineBiConsumer) {
        this.sampleName = Utils.nonNull(sampleName);
        this.records = Utils.nonNull(records).stream().sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR).collect(Collectors.toList());
        this.mandatoryColumns = Utils.nonNull(mandatoryColumns);
        this.dataLineToRecordFunction = Utils.nonNull(dataLineToRecordFunction);
        this.recordAndDataLineBiConsumer = Utils.nonNull(recordAndDataLineBiConsumer);
    }

    /**
     * Constructor given an input file, the mandatory column headers, and the lambdas for reading and writing records.
     * The sample name is read from a comment string in the file and the list of records is read using
     * the column headers and the appropriate lambda.
     * @throws UserException.BadInput if records are not sorted using {@link IntervalUtils#LEXICOGRAPHICAL_ORDER_COMPARATOR}
     */
    public TSVLocatableCollection(final File inputFile,
                                  final TableColumnCollection mandatoryColumns,
                                  final Function<DataLine, T> dataLineToRecordFunction,
                                  final BiConsumer<T, DataLine> recordAndDataLineBiConsumer) {
        IOUtils.canReadFile(inputFile);
        this.mandatoryColumns = Utils.nonNull(mandatoryColumns);
        this.dataLineToRecordFunction = Utils.nonNull(dataLineToRecordFunction);
        this.recordAndDataLineBiConsumer = Utils.nonNull(recordAndDataLineBiConsumer);

        try (final TSVReader reader = new TSVReader(inputFile)) {
            TableUtils.checkMandatoryColumns(reader.columns(), mandatoryColumns, UserException.BadInput::new);
            sampleName = reader.readSampleName();
            final List<T> recordsFromFile = reader.stream().collect(Collectors.toList());
            records = recordsFromFile.stream().sorted(IntervalUtils.LEXICOGRAPHICAL_ORDER_COMPARATOR).collect(Collectors.toList());
            if (!recordsFromFile.equals(records)) {
                throw new UserException.BadInput(String.format("Records in input file %s were not sorted in lexicographical order.", inputFile));
            }
        } catch (final IOException | UncheckedIOException e) {
            throw new UserException.CouldNotReadInputFile(inputFile, e);
        }
    }

    public String getSampleName() {
        return sampleName;
    }

    /**
     * @return  an unmodifiable view of the records contained in the collection
     */
    public List<T> getRecords() {
        return Collections.unmodifiableList(records);
    }

    /**
     * @return  a new modifiable list of {@link SimpleInterval}s corresponding to the {@link Locatable}s
     *          for each record contained in the collection
     */
    public List<SimpleInterval> getIntervals() {
        return records.stream()
                .map(r -> new SimpleInterval(r.getContig(), r.getStart(), r.getEnd()))
                .collect(Collectors.toList());
    }

    /**
     * Writes the sample name in a comment string and the records to file.
     */
    public void write(final File outputFile) {
        try (final TSVWriter writer = new TSVWriter(outputFile, sampleName)) {
            writer.writeSampleName();
            writer.writeAllRecords(records);
        } catch (final IOException e) {
            throw new UserException.CouldNotCreateOutputFile(outputFile, e);
        }
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) {
            return true;
        }
        if (o == null || getClass() != o.getClass()) {
            return false;
        }

        final TSVLocatableCollection<?> that = (TSVLocatableCollection<?>) o;

        return sampleName.equals(that.sampleName) && records.equals(that.records);
    }

    @Override
    public int hashCode() {
        int result = sampleName.hashCode();
        result = 31 * result + records.hashCode();
        return result;
    }

    private final class TSVReader extends TableReader<T> implements NamedSampleFile {
        private final File file;

        private TSVReader(final File file) throws IOException {
            super(file);
            this.file = file;
        }

        private String readSampleName() {
            return readSampleName(file);
        }

        @Override
        protected T createRecord(final DataLine dataLine) {
            Utils.nonNull(dataLine);
            return dataLineToRecordFunction.apply(dataLine);
        }
    }

    private final class TSVWriter extends TableWriter<T> {
        private final String sampleName;

        TSVWriter(final File file,
                  final String sampleName) throws IOException {
            super(file, mandatoryColumns);
            this.sampleName = Utils.nonNull(sampleName);
        }

        void writeSampleName() {
            try {
                writeComment(NamedSampleFile.SAMPLE_NAME_COMMENT_TAG + sampleName);
            } catch (final IOException e) {
                throw new UserException("Could not write sample name.");
            }
        }

        @Override
        protected void composeLine(final T record, final DataLine dataLine) {
            Utils.nonNull(record);
            Utils.nonNull(dataLine);
            recordAndDataLineBiConsumer.accept(record, dataLine);
        }
    }
}