package org.broadinstitute.hellbender.tools.copynumber.legacy.formats;

import org.broadinstitute.hellbender.utils.Utils;
import org.broadinstitute.hellbender.utils.mcmc.Decile;
import org.broadinstitute.hellbender.utils.mcmc.DecileCollection;
import org.broadinstitute.hellbender.utils.mcmc.ParameterEnum;
import org.broadinstitute.hellbender.utils.tsv.TableColumnCollection;

import java.io.File;
import java.util.AbstractMap;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.stream.Collectors;

/**
 * @author Samuel Lee &lt;slee@broadinstitute.org&gt;
 */
public final class ParameterDecileCollection<T extends Enum<T> & ParameterEnum> extends NamedSampleTSV<Map.Entry<T, DecileCollection>> {
    enum ParameterTableColumn {
        PARAMETER_NAME,
        POSTERIOR_10,
        POSTERIOR_20,
        POSTERIOR_30,
        POSTERIOR_40,
        POSTERIOR_50,
        POSTERIOR_60,
        POSTERIOR_70,
        POSTERIOR_80,
        POSTERIOR_90;

        static final TableColumnCollection COLUMNS = new TableColumnCollection((Object[]) values());
    }

    private final Map<T, DecileCollection> parameterToDecileCollectionMap;

    public ParameterDecileCollection(final String sampleName,
                                     final Map<T, DecileCollection> parameterToDecileCollectionMap,
                                     final Class<T> parameterClass,
                                     final String doubleFormat) {
        super(
                Utils.nonNull(sampleName),
                new ArrayList<>(parameterToDecileCollectionMap.entrySet()),
                ParameterTableColumn.COLUMNS,
                dataLine -> {
                    final String parameterName = dataLine.get(ParameterTableColumn.PARAMETER_NAME);
                    final T parameter = Enum.valueOf(Utils.nonNull(parameterClass), parameterName);
                    final DecileCollection deciles = new DecileCollection(Arrays.asList(
                            dataLine.getDouble(ParameterTableColumn.POSTERIOR_10),
                            dataLine.getDouble(ParameterTableColumn.POSTERIOR_20),
                            dataLine.getDouble(ParameterTableColumn.POSTERIOR_30),
                            dataLine.getDouble(ParameterTableColumn.POSTERIOR_40),
                            dataLine.getDouble(ParameterTableColumn.POSTERIOR_50),
                            dataLine.getDouble(ParameterTableColumn.POSTERIOR_60),
                            dataLine.getDouble(ParameterTableColumn.POSTERIOR_70),
                            dataLine.getDouble(ParameterTableColumn.POSTERIOR_80),
                            dataLine.getDouble(ParameterTableColumn.POSTERIOR_90)));
                    return new AbstractMap.SimpleEntry<>(parameter, deciles);},
                (record, dataLine) -> {
                    final T parameter = record.getKey();
                    final DecileCollection deciles = record.getValue();
                    dataLine.append(parameter.toString())
                            .append(String.format(doubleFormat, deciles.get(Decile.DECILE_10)))
                            .append(String.format(doubleFormat, deciles.get(Decile.DECILE_20)))
                            .append(String.format(doubleFormat, deciles.get(Decile.DECILE_30)))
                            .append(String.format(doubleFormat, deciles.get(Decile.DECILE_40)))
                            .append(String.format(doubleFormat, deciles.get(Decile.DECILE_50)))
                            .append(String.format(doubleFormat, deciles.get(Decile.DECILE_60)))
                            .append(String.format(doubleFormat, deciles.get(Decile.DECILE_70)))
                            .append(String.format(doubleFormat, deciles.get(Decile.DECILE_80)))
                            .append(String.format(doubleFormat, deciles.get(Decile.DECILE_90)));}
        );
        this.parameterToDecileCollectionMap = parameterToDecileCollectionMap;
    }

    public ParameterDecileCollection(final File file,
                                     final Class<T> parameterClass,
                                     final String doubleFormat) {
        super(
                Utils.nonNull(file),
                ParameterTableColumn.COLUMNS,
                dataLine -> {
                    final String parameterName = dataLine.get(ParameterTableColumn.PARAMETER_NAME);
                    final T parameter = Enum.valueOf(Utils.nonNull(parameterClass), parameterName);
                    final DecileCollection deciles = new DecileCollection(Arrays.asList(
                            dataLine.getDouble(ParameterTableColumn.POSTERIOR_10),
                            dataLine.getDouble(ParameterTableColumn.POSTERIOR_20),
                            dataLine.getDouble(ParameterTableColumn.POSTERIOR_30),
                            dataLine.getDouble(ParameterTableColumn.POSTERIOR_40),
                            dataLine.getDouble(ParameterTableColumn.POSTERIOR_50),
                            dataLine.getDouble(ParameterTableColumn.POSTERIOR_60),
                            dataLine.getDouble(ParameterTableColumn.POSTERIOR_70),
                            dataLine.getDouble(ParameterTableColumn.POSTERIOR_80),
                            dataLine.getDouble(ParameterTableColumn.POSTERIOR_90)));
                    return new AbstractMap.SimpleEntry<>(parameter, deciles);},
                (record, dataLine) -> {
                    final T parameter = record.getKey();
                    final DecileCollection deciles = record.getValue();
                    dataLine.append(parameter.toString())
                            .append(String.format(doubleFormat, deciles.get(Decile.DECILE_10)))
                            .append(String.format(doubleFormat, deciles.get(Decile.DECILE_20)))
                            .append(String.format(doubleFormat, deciles.get(Decile.DECILE_30)))
                            .append(String.format(doubleFormat, deciles.get(Decile.DECILE_40)))
                            .append(String.format(doubleFormat, deciles.get(Decile.DECILE_50)))
                            .append(String.format(doubleFormat, deciles.get(Decile.DECILE_60)))
                            .append(String.format(doubleFormat, deciles.get(Decile.DECILE_70)))
                            .append(String.format(doubleFormat, deciles.get(Decile.DECILE_80)))
                            .append(String.format(doubleFormat, deciles.get(Decile.DECILE_90)));}
        );
        parameterToDecileCollectionMap = getRecords().stream().collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue));
    }

    public DecileCollection getDeciles(final T parameter) {
        return parameterToDecileCollectionMap.get(parameter);
    }
}