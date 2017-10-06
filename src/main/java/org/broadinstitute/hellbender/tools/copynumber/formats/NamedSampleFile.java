package org.broadinstitute.hellbender.tools.copynumber.formats;

import org.broadinstitute.hellbender.exceptions.UserException;
import org.broadinstitute.hellbender.utils.io.IOUtils;
import org.broadinstitute.hellbender.utils.text.XReadLines;
import org.broadinstitute.hellbender.utils.tsv.TableUtils;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

public interface NamedSampleFile {
    String SAMPLE_NAME_COMMENT_TAG = "SAMPLE_NAME=";
    String SAMPLE_NAME_COMMENT_PREFIX = TableUtils.COMMENT_PREFIX + SAMPLE_NAME_COMMENT_TAG;

    default String readSampleName(final File file) {
        IOUtils.canReadFile(file);
        final List<String> sampleNameCommentLines = new ArrayList<>();
        try (final XReadLines reader = new XReadLines(file)) {
            for (final String line : reader) {
                if (!line.startsWith(SAMPLE_NAME_COMMENT_PREFIX)) {
                    break;
                }
                sampleNameCommentLines.add(line);
            }
        } catch (final IOException e) {
            throw new UserException.CouldNotReadInputFile(file);
        }
        if (sampleNameCommentLines.size() != 1) {
            throw new UserException.BadInput(String.format("File does not contain exactly one sample name specified by %s.",
                    SAMPLE_NAME_COMMENT_PREFIX));
        }
        return sampleNameCommentLines.get(0).replace(SAMPLE_NAME_COMMENT_PREFIX, "");
    }
}
